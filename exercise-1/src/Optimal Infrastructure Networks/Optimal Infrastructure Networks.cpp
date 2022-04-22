#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <math.h>

#define pair_index(id1, id2) id1 * num_cities + id2

const float DELTA = 0.1f;
const float GAMMA = 200.0f;

struct city_t {
	int id;
	std::string name;
	double population;
	float longitude;
	float latitude;
};

std::ostream& operator<<(std::ostream& os, city_t const& arg)
{
	os << arg.id << "(name = " << arg.name << ", pop = " << arg.population << ", lat = " << arg.latitude << ", lon = " << arg.longitude << ")";
	return os;
}

void floyd_warshall(std::vector<float> const& adj_matrix, std::vector<bool> const& util_matrix, std::vector<float>& mindist_matrix, size_t num_cities) {
	for (int i = 0; i < num_cities; i++) {
		mindist_matrix[pair_index(i, i)] = 0;
		for (int j = 0; j < num_cities; j++) {
			if (util_matrix[pair_index(i, j)]) {
				mindist_matrix[pair_index(i, j)] = adj_matrix[pair_index(i, j)];
			}
			else {
				mindist_matrix[pair_index(i, j)] = std::numeric_limits<float>::infinity();
			}
		}
	}

	for (int k = 0; k < num_cities; k++) {
		for (int i = 0; i < num_cities; i++) {
			if (i == k) continue;
			for (int j = i + 1; j < num_cities; j++) {
				if (j == k) continue;

				float route_over_k_length = mindist_matrix[pair_index(i, k)] + mindist_matrix[pair_index(k, j)];
				if (mindist_matrix[pair_index(i, j)] > route_over_k_length) {
					mindist_matrix[pair_index(i, j)] = route_over_k_length;
					mindist_matrix[pair_index(j, i)] = route_over_k_length;
				}
			}
		}
	}
}

float infra_cost(std::vector<float> const& adj_matrix, std::vector<bool> const& util_matrix, size_t num_cities) {
	float c = 0;

	for (int i = 0; i < num_cities; i++) {
		for (int j = i + 1; j < num_cities; j++) {
			int index = pair_index(i, j);
			if (util_matrix[index]) {
				c += adj_matrix[index];
			}
		}
	}

	return c;
}

float trans_cost(std::vector<float> const& mindist_matrix, std::vector<city_t> const& cities, size_t num_cities) {
	float c = 0;

	for (int i = 0; i < num_cities; i++) {
		for (int j = i + 1; j < num_cities; j++) {
			float effective_distance = (1.0f - DELTA) * mindist_matrix[pair_index(i, j)] + DELTA;
			c += cities[i].population * cities[j].population * effective_distance;
		}
	}

	return c; // / 2;
}

float get_total_cost(std::vector<float> const& adj_matrix, std::vector<bool> const& util_matrix, std::vector<float>& mindist_matrix, std::vector<city_t> const& cities, size_t num_cities) {
	floyd_warshall(adj_matrix, util_matrix, mindist_matrix, num_cities);

	float c_infra = infra_cost(adj_matrix, util_matrix, num_cities);
	float c_trans = trans_cost(mindist_matrix, cities, num_cities);

	return c_infra + GAMMA * c_trans;
}

std::vector<city_t> read_cities(std::ifstream &cities_file) {
	// read the city data
	std::vector<city_t> cities;

	std::string line;
	while (std::getline(cities_file, line)) {
		// id "name" pop lat lon
		size_t name_start = line.find("\t");
		size_t pop_start = line.find("\t", name_start + 1);
		size_t lat_start = line.find("\t", pop_start + 1);
		size_t lon_start = line.find("\t", lat_start + 1);

		int id = std::stoi(line.substr(0, name_start));
		std::string name = line.substr(name_start + 2, pop_start - 1 - (name_start + 2));
		double pop = std::stod(line.substr(pop_start + 1, lat_start - (pop_start + 1)));
		float lat = std::stof(line.substr(lat_start + 1, lon_start - (lat_start + 1)));
		float lon = std::stof(line.substr(lon_start + 1, line.length() - (lon_start + 1)));

		city_t city{ id, name, pop, lat, lon };

		cities.push_back(city);
	}

	return cities;
}

std::vector<float> read_adj_matrix(std::ifstream& distances_file, size_t num_cities) {
	// read the distances from the distances file
	std::vector<float> adj_matrix(num_cities * num_cities);

	std::string line;
	while (std::getline(distances_file, line)) {
		// id1 id2 dist
		size_t id2_start = line.find("\t");
		size_t dist_start = line.find("\t", id2_start + 1);

		int id1 = std::stoi(line.substr(0, id2_start));
		int id2 = std::stoi(line.substr(id2_start + 1, dist_start - (id2_start + 1)));
		float distance = std::stof(line.substr(dist_start + 1, line.length() - (dist_start + 1)));

		adj_matrix[pair_index(id1, id2)] = distance;
		adj_matrix[pair_index(id2, id1)] = distance;
	}

	return adj_matrix;
}

int main()
{
	// read cities
	std::ifstream cities_file;
	cities_file.open("cities.dat");
	if (!cities_file.is_open()) {
		std::cout << "Could not open cities.dat";
		return 1;
	}
	std::vector<city_t> cities = read_cities(cities_file);
	cities_file.close();

	size_t num_cities = cities.size();

	// read distances    
	std::ifstream distances_file;
	distances_file.open("distances.dat");
	if (!distances_file.is_open()) {
		std::cout << "Could not open distances_file.dat";
		return 1;
	}
	std::vector<float> adj_matrix = read_adj_matrix(distances_file, num_cities);
	distances_file.close();

	std::cout << "Finished reading " << num_cities << " cities." << std::endl;
 
	std::vector<bool> util_matrix(num_cities * num_cities, true);
	// TODO: init full, randomly, etc...
	std::vector<float> mindist_matrix(num_cities * num_cities);

	// random number generators
	std::random_device rd;
	std::mt19937 gen(rd());
	// for monte carlo rng
	std::uniform_real_distribution<> mc_dist(0, 1);
	// for edge selection
	std::uniform_int_distribution<> city_dist(0, num_cities - 1);

	float c_total_prev = get_total_cost(adj_matrix, util_matrix, mindist_matrix, cities, num_cities);;

	// number of iterations per temperature
	int its = 3 * num_cities * num_cities; 
	// initial temperature (steps: 10.0f, 1.0f, 0.1f, 0.01f)
	float T = 10.0f;

	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < its; i++) {
			// get random pair
			int id1 = city_dist(gen);
			int id2 = city_dist(gen);
			while (id1 == id2) {
				id2 = city_dist(gen);
			}

			int index1 = pair_index(id1, id2), index2 = pair_index(id2, id1);

			util_matrix[index1] = !util_matrix[index1];
			util_matrix[index2] = !util_matrix[index2];

			float c_total = get_total_cost(adj_matrix, util_matrix, mindist_matrix, cities, num_cities);
			if (c_total_prev > c_total || mc_dist(gen) > exp(-(c_total_prev - c_total) / T)) {
				c_total_prev = c_total;
			}
			else {
				util_matrix[index1] = !util_matrix[index1];
				util_matrix[index2] = !util_matrix[index2];
			}
		}

		T /= 10.0f;
	}


	// print results
	for (int i = 0; i < num_cities; i++) {
		for (int j = i + 1; j < num_cities; j++) {
			if (util_matrix[pair_index(i, j)]) {
				std::cout << cities[i].latitude << " " << cities[i].longitude << " " << cities[j].latitude << " " << cities[j].longitude << std::endl;
			}
		}
	}

	std::cout << "Done" << std::endl;
}