/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  //Set the number of particles.
  this->num_particles = 10000;
  this->particles.resize(this->num_particles);
  this->weights.resize(this->num_particles);

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
    this->weights[i] = 1.0;
  }

  this->is_initialized = true;
#if PRINT
  printParticles();
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
  std::cout << "Prediction step\n";
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
  printParticles();
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {

  std::cout << "Update step\n";
	//Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

  double gaussFactor = 1.0/(2*M_PI*std_landmark[0]*std_landmark[1]);
  double sig_x_sqr = std_landmark[0]*std_landmark[0];
  double sig_y_sqr = std_landmark[1]*std_landmark[1];

  for(size_t i=0; i < this->num_particles; ++i)
  {
    double x = this->particles[i].x;
    double y = this->particles[i].y;
    double theta = this->particles[i].theta;

    //Transform observations in map space
    auto transformed = observations;
    for(size_t j=0; j < transformed.size(); ++j)
    {
      //Translation
      double x_trans = x + transformed[j].x;
      double y_trans = y + transformed[j].y;

      //Rotation
      double cos_theta = std::cos(theta);
      double sin_theta = std::sin(theta);
      transformed[j].x = x_trans*cos_theta - y_trans*sin_theta;
      transformed[j].y = x_trans*sin_theta + y_trans*cos_theta;

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

//      std::cout << "Distance: x: " << transformed[j].x - map_landmarks.landmark_list[transformed[j].id-1].x_f
//          << " y: " << transformed[j].y - map_landmarks.landmark_list[transformed[j].id-1].y_f << std::endl;
      //Update weights.
      double x_val = std::pow(transformed[j].x - map_landmarks.landmark_list[transformed[j].id-1].x_f, 2)/(2*sig_x_sqr);
      double y_val = std::pow(transformed[j].y - map_landmarks.landmark_list[transformed[j].id-1].y_f, 2)/(2*sig_y_sqr);
      double weight = gaussFactor*std::exp(-1.0*(x_val+y_val));
      //std::cout << "Weight update " << weight << std::endl;
      if(weight > 0.0)
        this->particles[i].weight *= weight;
    }

  }

#if PRINT
  printParticles();
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
  std::cout << "Resample step\n";
  std::default_random_engine generator;
  std::uniform_int_distribution<size_t> uni_dist(0, this->num_particles-1);

  size_t idx = uni_dist(generator);

  auto resampledParticles = this->particles;

  //Find max weigth
  double maxWeight = 0;
  for(size_t i=0; i < this->num_particles; ++i)
  {
    if(particles[i].weight > maxWeight)
      maxWeight = particles[i].weight;
  }

  std::cout << "MAX weight " << maxWeight << std::endl;
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

  //  std::cout << "idx " << idx << std::endl;
    resampledParticles[i] = particles[idx%this->num_particles];

  }

  particles = resampledParticles;

#if PRINT
  printParticles();
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
void ParticleFilter::printParticles()
{
  for(size_t i=0; i < this->num_particles; ++i)
  {
    std::cout << "[" << this->particles[i].id << "]"
        << " x: " << this->particles[i].x
        << " y: " << this->particles[i].y
        << " w: " << this->particles[i].weight
        << std::endl;
  }
}
