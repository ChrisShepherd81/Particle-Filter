/*
 * particle_filter.cpp
 *
 *  Created on: May 02, 2017
 *      Author: christian@inf-schaefer.de
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

#define PRINT 0
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::init(double x, double y, double theta, double std[])
{
  //Set the number of particles.
  this->num_particles = 100;
  this->particles.resize(this->num_particles);

  std::default_random_engine generator;
  std::normal_distribution<double> dist_x(x,std[0]);
  std::normal_distribution<double> dist_y(y,std[1]);
  std::normal_distribution<double> dist_theta(theta,std[2]);

  //Initialize all particles to first position (based on estimates of
  //x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  for(size_t i=0; i < this->num_particles; i++)
  {
    this->particles[i].id = i;
    this->particles[i].x = dist_x(generator);
    this->particles[i].y = dist_y(generator);
    this->particles[i].theta = dist_theta(generator);
    this->particles[i].weight = 1.0;
  }

  this->is_initialized = true;

#if PRINT
  printParticles("Init");
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
  std::default_random_engine generator;
  std::normal_distribution<double> dist_x(0.0,std_pos[0]);
  std::normal_distribution<double> dist_y(0.0,std_pos[1]);
  std::normal_distribution<double> dist_theta(0.0,std_pos[2]);

  double y_diff_yaw_rate, d_theta;
  if(std::abs(yaw_rate) > 1e-6)
  {
    y_diff_yaw_rate = velocity/yaw_rate;
    d_theta = delta_t*yaw_rate;
  }
  else
  {
    d_theta = 0.0;
  }

  for(size_t i=0; i < this->num_particles; ++i)
  {
    double theta = this->particles[i].theta;
    double new_x, new_y;

    if(std::abs(yaw_rate) > 1e-6)
    {
      new_x = this->particles[i].x + (y_diff_yaw_rate*((std::sin(theta+d_theta)-std::sin(theta))));
      new_y = this->particles[i].y + (y_diff_yaw_rate*((std::cos(theta)-std::cos(theta+d_theta))));
    }
    else
    {
      new_x = this->particles[i].x + (velocity*(std::cos(theta)));
      new_y = this->particles[i].y + (velocity*(std::sin(theta)));
    }

    double new_theta = theta+d_theta;

    // Add measurements to each particle and add random Gaussian noise.
    this->particles[i].x = new_x + dist_x(generator);
    this->particles[i].y = new_y + dist_y(generator);
    this->particles[i].theta = new_theta + dist_theta(generator);
  }

#if PRINT
  printParticles("Prediction");
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks)
{
	//Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

  double gaussFactor = 1.0/(2*M_PI*std_landmark[0]*std_landmark[1]);
  double sig_x_sqr = std_landmark[0]*std_landmark[0];
  double sig_y_sqr = std_landmark[1]*std_landmark[1];

#if 1 //PRINT
  std::cout << "Gauss Factor: " << gaussFactor << std::endl;
  std::cout << "(Sigma x)^2: " << sig_x_sqr << std::endl;
  std::cout << "(Sigma y)^2: " << sig_y_sqr << std::endl;
#endif

  static size_t timestamp = 1;
  std::cout << "Timestamp " << timestamp++ << std::endl;

  for (size_t j=0; j < observations.size(); ++j)
  {
    std::cout << "Landmark " << j << "(" << observations[j].x << "," << observations[j].y << ")\n";
  }

  for(size_t i=0; i < this->num_particles; ++i)
  {
    double p_x = this->particles[i].x;
    double p_y = this->particles[i].y;
    double theta = this->particles[i].theta;

    //Transform observations in map space
    auto transformed = observations;
    for(size_t j=0; j < transformed.size(); ++j)
    {
      //Rotation and translation
      double cos_theta = std::cos(theta);
      double sin_theta = std::sin(theta);
      transformed[j].x = observations[j].x*cos_theta - observations[j].y*sin_theta + p_x;
      transformed[j].y = observations[j].x*sin_theta + observations[j].y*cos_theta + p_y;

      //Associate transformed observations to landmarks with nearest neighbor
      double min = sensor_range;
      for(size_t k=0; k < map_landmarks.landmark_list.size(); ++k)
      {
        double d = dist(transformed[j].x, transformed[j].y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
        if(d < std::abs(min))
        {
          min = d;
          transformed[j].id = map_landmarks.landmark_list[k].id_i;
        }
      }

      //Update weights.
      double x_val = std::pow(transformed[j].x - map_landmarks.landmark_list[transformed[j].id-1].x_f, 2)/(2*sig_x_sqr);
      double y_val = std::pow(transformed[j].y - map_landmarks.landmark_list[transformed[j].id-1].y_f, 2)/(2*sig_y_sqr);
      double weight = gaussFactor*std::exp(-(x_val+y_val));

      if(weight > 0.0)
      {
        std::cout << "Weight is "  << weight << std::endl;
        this->particles[i].weight *= weight;
        if(isnan(this->particles[i].weight))
        {
          std::cout << "FATAL weight is "  << weight;
          exit(1);
        }
      }
    }
  }

#if 1 //PRINT
  printParticles("Update");
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::resample()
{
	//Resample particles with replacement with probability proportional to their weight.
  std::default_random_engine generator;
  std::uniform_int_distribution<size_t> uni_dist(0, this->num_particles-1);

  size_t idx = uni_dist(generator);

  auto resampledParticles = this->particles;

  //Calculate weight sum
  double weight_sum = 0;
  for(size_t i=0; i < this->num_particles; ++i)
  {
   	weight_sum += particles[i].weight;
  }

  //Normalize weights
  for(size_t i=0; i < this->num_particles; ++i)
  {
   	particles[i].weight = particles[i].weight / weight_sum;
  }

  //Find max weigth
  double maxWeight = 0;
  for(size_t i=0; i < this->num_particles; ++i)
  {
    if(particles[i].weight > maxWeight)
      maxWeight = particles[i].weight;
  }

  std::uniform_real_distribution<double> real_uni_dist(0, 2*maxWeight);
  double beta = 0;

  for(size_t i=0; i < this->num_particles; ++i)
  {
    beta += real_uni_dist(generator);

    while(particles[idx].weight < beta)
    {
      beta -= particles[idx].weight;
      idx = (idx+1)%this->num_particles ;
    }

    resampledParticles[i] = particles[idx];
    resampledParticles[i].id = i+1;
  }

  particles = resampledParticles;
  write("../data/particles.txt");

#if 1 //PRINT
  printParticles("Resample");
#endif

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::printParticles(std::string action)
{
  std::cout << action << ":\n";
  for(size_t i=0; i < this->num_particles; ++i)
  {
    std::cout << "[" << this->particles[i].id << "]"
        << " x: " << this->particles[i].x
        << " y: " << this->particles[i].y
        << " theta: " << this->particles[i].theta
        << " w: " << this->particles[i].weight
        << std::endl;
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  //Clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

  return particle;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string ParticleFilter::getSenseX(Particle best)
{
  vector<double> v = best.sense_x;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string ParticleFilter::getSenseY(Particle best)
{
  vector<double> v = best.sense_y;
  stringstream ss;
  copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
