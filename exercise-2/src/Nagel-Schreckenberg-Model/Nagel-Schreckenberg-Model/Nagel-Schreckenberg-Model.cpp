#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <math.h>
#include <format>

void update_nagel_schreckenberg(std::vector<int>& positions, std::vector<int>& velocities, int length, int vmax, float p, std::uniform_real_distribution<float> &dec_distribution, std::mt19937 & gen)
{
	for (size_t i = 0; i < positions.size(); i++)
	{
		// Acceleration: vi → vi + 1 if vi < vmax
		if (velocities[i] < vmax) velocities[i]++;
		//else {
		//	float thr = exp(-(velocities[i] - 4.0f));
		//	if (dec_distribution(gen) < thr) {
		//		velocities[i]++;
		//	}
		//}

		// Collision avoidance : vi → max(vi, x(i + 1) − xi − 1)
		int next_car_i = (i + 1) % positions.size();
		int collision_velocity = (positions[next_car_i] - positions[i] - 1 + length) % length;
		if (velocities[i] > collision_velocity) velocities[i] = collision_velocity;

		// Random deceleration : vi → vi − 1 with probability p if vi > 0
		float threshold = dec_distribution(gen);
		if (velocities[i] > 0 && threshold < p) velocities[i]--;
	}

	// Movement : xi → xi + vi
	for (size_t i = 0; i < positions.size(); i++)
	{
		positions[i] = (positions[i] + velocities[i]) % length;
	}
}

std::vector<int> prepare_positions(int cars, int length)
{
	std::vector<int> positions(cars);

	std::vector<int> possible_positions(length);
	for (int l = 0; l < length; l++)
	{
		possible_positions[l] = l;
	}

	for (int i = 0; i < cars; i++)
	{
		int position_index = rand() % possible_positions.size();

		positions[i] = possible_positions[position_index];
		std::swap(possible_positions[position_index], possible_positions.back());
		possible_positions.pop_back();
	}

	std::sort(positions.begin(), positions.end());

	return positions;
}

int main()
{
	int length = 1000;
	int vmax = 5;
	float p = 0.35f;
	
	int time_steps = 1000;

	float density_min = 0.0f;
	float density_max = 0.6f;
	int density_steps = 100;

	// random number generators
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dec_distribution(0.0f, 1.0f);


	for (int d = 0; d < density_steps; d++)
	{
		float desired_density = density_min + (density_max - density_min) / (density_steps - 1) * d;
		int cars = (int)(length * desired_density);
		float real_density = (float)cars / length;

		std::cout << "Computing density rho=" << real_density;

		std::vector<int> positions = prepare_positions(cars, length);
		std::vector<int> velocities(cars, 0);

		std::string res = "T;P;V\n";

		double flow = 0;

		for (int t = 0; t < time_steps; t++)
		{
			for (int i = 0; i < cars; i++)
			{
				res += std::to_string(t) + ";" + std::to_string(positions[i]) + ";" + std::to_string(velocities[i]) + "\n";
				flow += velocities[i];
			}

			update_nagel_schreckenberg(positions, velocities, length, vmax, p, dec_distribution, gen);
		}

		std::cout << " flow=" << flow / time_steps / length;

		std::string filename = std::format("tpv{}.csv", d);
		std::cout << " writing to " << filename;

		std::ofstream f;
		f.open(filename);
		f << res;
		f.close();

		std::cout << " done" << std::endl;
	}

	return 0;

}
