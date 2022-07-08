#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>    
#include <stdexcept>
#include <numeric>
#include <string>
#include <execution>

struct shared_ride_distances {
	bool focal_user_was_paired;
	float focal_user_detour;
	float total_distance;
};

float *distance_lut;

void generate_distance_lut(int M) {
	distance_lut = new float[2 * M + 1];
	for (int m = -M; m <= M; m++) {
		distance_lut[m + M] = (float)std::sqrt(2 - 2 * std::cos(2 * M_PI * m / M));
	}
}

inline float get_distance(int M, int m1, int m2)
{
	return distance_lut[m1 - m2 + M];
}

inline float get_utility_diff(float beta, float detour)
{
	return 1 - beta * detour;
}

void update_adoption_probabilities(float dt, std::vector<float> & p, const std::vector<float> & expected_utility_diffs)
{
	for (int i = 0; i < p.size(); i++) 
	{
		p[i] = p[i] + dt * p[i] * (1 - p[i]) * expected_utility_diffs[i];
		// TODO: Wertebereich einschränken
		p[i] = std::min(std::max(p[i], 1e-3f), 1.0f - 1e-3f);
	}
}

int get_offset_skipped_index(int i, int requests, int offset, int skip)
{
	if (i >= skip) {
		offset++;
	}
	return (i + offset) % requests;
}

shared_ride_distances compute_shared_ride_distances(
	int M, int requests, int offset, int skip,
	std::vector<int>& destination_buffer, std::vector<int>& idx_buffer) {

	// no need to consider solo rides, since they always contribute the same constant distance 2
	float total_distance = 0;
	float focal_user_detour = -1.0f;
	bool focal_user_was_paired = false;


	// iterate only the pairs
	int pairs = requests / 2;
	for (int i = 0; i < pairs; i++)
	{
		int i1 = get_offset_skipped_index(2 * i + 0, requests, offset, skip);
		int i2 = get_offset_skipped_index(2 * i + 1, requests, offset, skip);
		float d = get_distance(M, destination_buffer[idx_buffer[i1]], destination_buffer[idx_buffer[i2]]);

		if (idx_buffer[i1] == 0 || idx_buffer[i2] == 0)
		{
			focal_user_was_paired = true;
			focal_user_detour = d;
		}

		total_distance += d;
	}

	return { focal_user_was_paired, focal_user_detour, total_distance };
}

void update_shortest_shared_ride(
	shared_ride_distances& shortest_ride,
	int M, int requests, int offset, int skip,
	std::vector<int>& destination_buffer, std::vector<int>& idx_buffer)
{
	shared_ride_distances d = compute_shared_ride_distances(M, requests, offset, skip, destination_buffer, idx_buffer);
	if (d.total_distance < shortest_ride.total_distance)
	{
		shortest_ride = d;
	}
}

shared_ride_distances get_shortest_shared_ride(
	int M, int requests,
	std::vector<int>& destination_buffer, std::vector<int>& idx_buffer) {
	if (requests < 2) {
		return { true, 0, 0 };
	}

	// no need to consider solo rides, since they always contribute the same constant distance 2
	shared_ride_distances shortest_ride = { false, -1.0f, 1e10 };

	for (int offset = 0; offset < 2; offset++) {
		if (requests % 2 == 0)
		{
			update_shortest_shared_ride(shortest_ride, M, requests, offset, 1000, destination_buffer, idx_buffer);
		}
		else
		{
			for (int skip = 0; skip < requests; skip++)
			{
				update_shortest_shared_ride(shortest_ride, M, requests, offset, skip, destination_buffer, idx_buffer);
			}
		}
	}

	return shortest_ride;
}

inline int get_random_destination(int M, std::function<int()> rng)
{
	return rng();
}

void update_approximated_utility_diffs(
	int realizations, float beta, int N, int M, float dt, const std::vector<float> & p,
	std::vector<int>& destination_buffer, std::vector<int>& idx_buffer,
	std::vector<float>& utility_diff_buffer, std::vector<int>& destination_counter_buffer,
	std::function<int()> rng, std::function<float()> rng_float)
{
	for (int i = 0; i < M; i++)
	{
		utility_diff_buffer[i] = 0;
		destination_counter_buffer[i] = 0;
	}

	for (int i = 0; i < realizations; i++)
	{
		// compute how many people actually share a ride
		int ride_sharing_requests = 1;
		destination_buffer[0] = get_random_destination(M, rng);

		for (int n = 1; n < N; n++) {
			int destination = get_random_destination(M, rng);
			if (rng_float() < p[destination]) {
				destination_buffer[ride_sharing_requests] = destination;
				ride_sharing_requests++;
			}
		}

		// select random destinations, first one is fixed
		for (int n = 0; n < ride_sharing_requests; n++)
		{
			idx_buffer[n] = n;
		}

		// sort indecies
		std::stable_sort(idx_buffer.begin(), idx_buffer.begin() + ride_sharing_requests,
			[&destination_buffer](size_t i1, size_t i2) {return destination_buffer[i1] < destination_buffer[i2]; });

		// get shortest shared ride
		shared_ride_distances shortest_shared_ride = get_shortest_shared_ride(M, ride_sharing_requests, destination_buffer, idx_buffer);

		if (shortest_shared_ride.focal_user_was_paired) {
			int focal_destination_count = 1;
			for (int n = 1; n < ride_sharing_requests; n++)
			{
				if (destination_buffer[n] == destination_buffer[0])
				{
					focal_destination_count++;
				}
			}

			utility_diff_buffer[destination_buffer[0]] += get_utility_diff(beta, 0.5f * shortest_shared_ride.focal_user_detour / focal_destination_count);
			destination_counter_buffer[destination_buffer[0]]++;
		}
		else {
			// was not assigned a partner and therefore took the trip alone
			utility_diff_buffer[destination_buffer[0]] += 0;
		}
	}

	for (int i = 0; i < M; i++)
	{
		utility_diff_buffer[i] /= destination_counter_buffer[i];
	}
}

struct ride_sharing_time_series {
	std::vector<float> acc_probabilites_time_series;
	std::vector<float> destination_probabilities;
};

ride_sharing_time_series get_ride_sharing_time_series(
	int realizations, int time_steps, float dt,
	float beta, int M, int N, float p0,
	std::function<int()> rng,
	std::function<float()> rng_float
)
{
	std::vector<int> destination_buffer = std::vector<int>(N);
	std::vector<int> idx_buffer = std::vector<int>(N);

	std::vector<float> utility_diff_buffer = std::vector<float>(M);
	std::vector<int> destination_counter_buffer = std::vector<int>(M);

	std::vector<float> acc_probabilites_time_series = std::vector<float>(time_steps);
	std::vector<float> p = std::vector<float>(M, p0);

	for (int s = 0; s < time_steps; s++)
	{
		update_approximated_utility_diffs(
			realizations, beta, N, M, dt, p, destination_buffer, idx_buffer, utility_diff_buffer, destination_counter_buffer, rng, rng_float
		);
		update_adoption_probabilities(dt, p, utility_diff_buffer);

		acc_probabilites_time_series[s] = std::accumulate(p.begin(), p.end(), 0.0f) / M;
	}

	return { acc_probabilites_time_series, p };
}

void save_ride_sharing_time_series(ride_sharing_time_series ts, int N, int id)
{
	std::string fileid = std::to_string(N) + "_" + std::to_string(id) + ".csv";
	std::cout << " writing to *" << fileid << " ... ";

	std::ofstream f_ts;
	f_ts.open("ts_" + fileid);
	f_ts << "step;p" << std::endl;
	for (int i = 0; i < ts.acc_probabilites_time_series.size(); i++)
	{
		f_ts << i << ";" << ts.acc_probabilites_time_series[i] << std::endl;
	}
	f_ts.close();

	std::ofstream f_dp;
	f_dp.open("dp_" + fileid);
	f_dp << "m;p" << std::endl;
	for (int i = 0; i < ts.destination_probabilities.size(); i++)
	{
		f_dp << i << ";" << ts.destination_probabilities[i] << std::endl;
	}
	f_dp.close();

	std::cout << "done" << std::endl;
}

int main()
{
	std::cout << "Hello World!" << std::endl;

	int realizations = 1000;
	int time_steps = 1000;
	float dt = 0.1f;
	int M = 360;
	int N = 2;

	float p0 = 0.5f;

	float beta_min = 1.0f;
	float beta_max = 3.0f;
	int beta_steps = 21;

	generate_distance_lut(M);

	std::vector<int> b = std::vector<int>(beta_steps);
	std::iota(b.begin(), b.end(), 0);

	std::for_each(
		std::execution::par,
		b.begin(),
		b.end(),
		[&beta_min, &beta_max, &beta_steps, &time_steps, &dt, &N, &M, &p0, &realizations](int item) {
			float beta = beta_min + (beta_max - beta_min) * item / (beta_steps - 1);
			std::stringstream msg;
			msg << item << ": computing N=" << N << ", p0=" << p0 << "(realizations=" << realizations << ", time_steps=" << time_steps << ") at beta=" << beta << std::endl;
			std::cout << msg.str();

			// random number generation
			std::random_device dev;
			std::mt19937 rng(dev());
			std::uniform_int_distribution<> dist(0, M - 1);
			std::uniform_real_distribution<> dist_float(0, 1);

			ride_sharing_time_series ts = get_ride_sharing_time_series(
				realizations, time_steps, dt, beta, M, N, p0,
				[&dist, &rng]() -> int {
					return dist(rng);
				},
				[&dist_float, &rng]() -> float {
					return (float)dist_float(rng);
				}
				);
			save_ride_sharing_time_series(ts, N, item);
		}
	);
}
