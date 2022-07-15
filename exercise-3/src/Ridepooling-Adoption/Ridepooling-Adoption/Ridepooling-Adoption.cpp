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
	float focal_user_detour;
	float total_distance;
};

inline float get_distance(int M, int m1, int m2, const std::vector<float>& distance_lut)
{
	return distance_lut[m1 - m2 + M];
}

inline float get_utility_diff(float beta, float detour)
{
	return 1 - beta * detour;
}

inline float clamp(float v, float min, float max)
{
	if (v < min) return min;
	if (v > max) return max;
	return v;

	// alternative: return std::min(std::max(v, min), max);
}

void update_adoption_probabilities(float dt, std::vector<float> & p, const std::vector<float> & expected_utility_diffs)
{
	for (int i = 0; i < p.size(); i++) 
	{
		// Wertebereich einschränken
		p[i] = clamp(p[i] + dt * p[i] * (1 - p[i]) * expected_utility_diffs[i], 1e-3f, 1.0f - 1e-3f);
	}
}

int get_offset_skipped_index(int i, int requests, int offset, int skip)
{
	if (i >= skip) {
		offset++;
	}
	return (i + offset) % requests;
}

int c = 0;

void perfect_matching_branch_and_bound(
	int M, int requests, int current, float sum, float focal_user_detour,
	const std::vector<float>& distance_lut, std::vector<bool>& used_buffer,
	const std::vector<int>& destination_buffer, const std::vector<int>& idx_buffer,
	shared_ride_distances& shortest_shared_ride) {

	if (sum >= shortest_shared_ride.total_distance) return;

	if (current == requests) {
		c++;
		shortest_shared_ride = { focal_user_detour, sum };
		return;
	}

	if (used_buffer[current])
	{
		perfect_matching_branch_and_bound(M, requests, current + 1, sum, focal_user_detour, distance_lut, used_buffer, destination_buffer, idx_buffer, shortest_shared_ride);
	}
	else
	{
		used_buffer[current] = true;
		for (int i = 0; i < requests; i++)
		{
			if (!used_buffer[i])
			{
				used_buffer[i] = true;

				int m1 = destination_buffer[idx_buffer[current]];
				int m2 = destination_buffer[idx_buffer[i]];

				float d = get_distance(M, m1, m2, distance_lut);
				if (idx_buffer[current] == 0) {
					// is focal user
					focal_user_detour = d;
				}
				perfect_matching_branch_and_bound(M, requests, current + 1, sum + d, focal_user_detour, distance_lut, used_buffer, destination_buffer, idx_buffer, shortest_shared_ride);
				used_buffer[i] = false;
			}
		}
		used_buffer[current] = false;
	}
}

shared_ride_distances compute_shared_ride_distances(
	int M, int requests, int offset, int skip, const std::vector<float>& distance_lut,
	const std::vector<int>& destination_buffer, const std::vector<int>& idx_buffer
) {
	float total_distance = 0.0f;
	float focal_user_detour = 0.0f;

	// iterate only the pairs
	int pairs = requests / 2;
	for (int i = 0; i < pairs; i++)
	{
		int i1 = get_offset_skipped_index(2 * i + 0, requests, offset, skip);
		int i2 = get_offset_skipped_index(2 * i + 1, requests, offset, skip);
		float d = get_distance(M, destination_buffer[idx_buffer[i1]], destination_buffer[idx_buffer[i2]], distance_lut);

		if (idx_buffer[i1] == 0 || idx_buffer[i2] == 0)
		{
			focal_user_detour = d;
		}

		total_distance += d;
	}

	return { focal_user_detour, total_distance };
}

void update_shortest_shared_ride(
	shared_ride_distances& shortest_ride, const std::vector<float>& distance_lut,
	int M, int requests, int offset, int skip,
	const std::vector<int>& destination_buffer, const std::vector<int>& idx_buffer)
{
	shared_ride_distances d = compute_shared_ride_distances(M, requests, offset, skip, distance_lut, destination_buffer, idx_buffer);
	if (d.total_distance < shortest_ride.total_distance)
	{
		shortest_ride = d;
	}
}

inline shared_ride_distances get_shortest_shared_ride_general(
	int N, int M, int requests, const std::vector<float>& distance_lut,
	std::vector<int> & destination_buffer, std::vector<int> & idx_buffer
) {
	if (requests < 2) {
		return { 0, 0 };
	}

	if (N == 2) 
	{
		float d = get_distance(M, destination_buffer[0], destination_buffer[1], distance_lut);
		return { d, d };
	}

	std::sort(idx_buffer.begin(), idx_buffer.begin() + requests,
		[&destination_buffer](size_t i1, size_t i2) { return destination_buffer[i1] < destination_buffer[i2]; }
	);

	// no need to consider solo rides, since they always contribute the same constant distance 2
	shared_ride_distances shortest_ride = { -1.0f,  std::numeric_limits<float>::infinity() };

	for (int offset = 0; offset < 2; offset++) {
		if (requests % 2 == 0)
		{
			update_shortest_shared_ride(shortest_ride, distance_lut, M, requests, offset, 1000,  destination_buffer, idx_buffer);
		}
		else
		{
			for (int skip = 0; skip < requests; skip++)
			{
				update_shortest_shared_ride(shortest_ride, distance_lut, M, requests, offset, skip, destination_buffer, idx_buffer);
			}
		}
	}

	return shortest_ride;
}

shared_ride_distances perfect_matching_branch_and_bound(
	int N, int M, int requests,
	const std::vector<float>& distance_lut, std::vector<bool>& used_buffer,
	std::vector<int>& destination_buffer, std::vector<int>& idx_buffer) {

	for (int i = 0; i < requests; i++)
	{
		used_buffer[i] = false;
	}

	shared_ride_distances estimate = get_shortest_shared_ride_general(N, M, requests, distance_lut, destination_buffer, idx_buffer);

	perfect_matching_branch_and_bound(M, requests, 0, 0.0f, 0.0f, distance_lut, used_buffer, destination_buffer, idx_buffer, estimate);

	return estimate;
}

inline int get_random_destination(std::mt19937& rng, std::uniform_int_distribution<int>& dist_int_0_Mm1)
{
	return dist_int_0_Mm1(rng);
}

void update_approximated_utility_diffs(
	int realizations, float beta, int N, int M, float dt, const std::vector<float> & p, const std::vector<float> &distance_lut,
	std::vector<int>& destination_buffer, std::vector<int>& idx_buffer, std::vector<bool> & used_buffer,
	std::vector<float>& utility_diff_buffer, std::vector<int>& destination_counter_buffer,
	std::mt19937 & rng, 
	std::uniform_real_distribution<float> & dist_float_0_1, std::uniform_int_distribution<int> & dist_int_0_Mm1)
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
		int focal_destination_count = 1;
		destination_buffer[0] = get_random_destination(rng, dist_int_0_Mm1);

		for (int n = 1; n < N; n++) {
			int destination = get_random_destination(rng, dist_int_0_Mm1);
			if (dist_float_0_1(rng) < p[destination]) {
				destination_buffer[ride_sharing_requests] = destination;
				if (destination == destination_buffer[0])
				{
					focal_destination_count++;
				}

				ride_sharing_requests++;
			}
		}

		// select random destinations, first one is fixed
		for (int n = 0; n < ride_sharing_requests; n++)
		{
			idx_buffer[n] = n;
		}

		// get shortest shared ride
		shared_ride_distances shortest_shared_ride = perfect_matching_branch_and_bound(N, M, ride_sharing_requests, distance_lut, used_buffer, destination_buffer, idx_buffer);

		utility_diff_buffer[destination_buffer[0]] += get_utility_diff(beta, 0.5f * shortest_shared_ride.focal_user_detour / focal_destination_count);
		destination_counter_buffer[destination_buffer[0]]++;
	}

	for (int i = 0; i < M; i++)
	{
		if (destination_counter_buffer[i] != 0)
			utility_diff_buffer[i] /= destination_counter_buffer[i];
	}
}

struct ride_sharing_time_series {
	std::vector<float> acc_probabilites_time_series;
	std::vector<float> destination_probabilities;
};

void save_probability_snapshot(std::vector<float>& p, int N, int id, int s)
{
	//std::string fileid = std::to_string(N) + "_" + std::to_string(id) + "-" + std::to_string(s) + ".csv";
	//std::cout << " writing to *" << fileid << " ... ";

	//std::ofstream f_dp;
	//f_dp.open("dp_" + fileid);
	//f_dp << "m;p" << std::endl;
	//for (int i = 0; i < p.size(); i++)
	//{
	//	f_dp << i << ";" << p[i] << std::endl;
	//}
	//f_dp.close();

	//std::cout << "done" << std::endl;
	std::cout << s << std::endl;
}

ride_sharing_time_series get_ride_sharing_time_series(
	int realizations, int time_steps, float dt,
	float beta, int M, int N, float p0, int seed)
{
	// random number generation
	std::mt19937 rng(seed);
	std::uniform_real_distribution<float> dist_float_0_1(0, 1);
	std::uniform_int_distribution<int> dist_int_0_Mm1(0, M - 1);

	std::vector<int> destination_buffer = std::vector<int>(N);
	std::vector<int> idx_buffer = std::vector<int>(N);
	std::vector<bool> used_buffer = std::vector<bool>(N);

	std::vector<float> utility_diff_buffer = std::vector<float>(M);
	std::vector<int> destination_counter_buffer = std::vector<int>(M);

	std::vector<float> acc_probabilites_time_series = std::vector<float>(time_steps);
	std::vector<float> p = std::vector<float>(M, p0);

	// generate distance lut
	std::vector<float> distance_lut = std::vector<float>(2 * M + 1);
	for (int m = -M; m <= M; m++) {
		distance_lut[m + M] = (float)std::sqrt(2 - 2 * std::cos(2 * M_PI * m / M));
	}

	for (int s = 0; s < time_steps; s++)
	{
		if (s % 100 == 0) save_probability_snapshot(p, N, 5, s / 50);
		update_approximated_utility_diffs(
			realizations, beta, N, M, dt, p, distance_lut, destination_buffer, idx_buffer, used_buffer, utility_diff_buffer, destination_counter_buffer, rng, dist_float_0_1, dist_int_0_Mm1
		);
		update_adoption_probabilities(dt, p, utility_diff_buffer);

		acc_probabilites_time_series[s] = std::accumulate(p.begin(), p.end(), 0.0f) / M;
		if (isnan(acc_probabilites_time_series[s]))
		{
			std::cout << "error";
		}
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

	int realizations = 10000;
	int time_steps = 10000;
	float dt = 0.1f;
	int M = 360;
	int N = 8;

	float p0 = 0.05f;

	float beta_min = 2.0f;
	float beta_max = 5.0f;
	int beta_steps = 31;

	std::vector<int> b = std::vector<int>(beta_steps);
	std::iota(b.begin(), b.end(), 0);

	std::for_each(
		std::execution::par,
		b.begin(),
		b.end(),
		[&beta_min, &beta_max, &beta_steps, &time_steps, &dt, &N, &M, &p0, &realizations](int item) {
			float beta = beta_min +(beta_max - beta_min) * item / (beta_steps - 1);
			std::stringstream msg;
			msg << item << ": computing N=" << N << ", p0=" << p0 << "(realizations=" << realizations << ", time_steps=" << time_steps << ") at beta=" << beta << std::endl;
			std::cout << msg.str();

			int seed = time(nullptr);

			ride_sharing_time_series ts = get_ride_sharing_time_series(
				realizations, time_steps, dt, beta, M, N, p0, seed
			);
			save_ride_sharing_time_series(ts, N, item + 16);
		}
	);

	std::cout << "IN " << c << " OF " << time_steps * realizations << " USED";
}
