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
  this->num_particles = 100;
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

  for(size_t i=0; i < this->num_particles; i++)
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
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
