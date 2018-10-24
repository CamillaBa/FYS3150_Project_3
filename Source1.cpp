#include <iostream>
#include <vector> 
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip> 

double ** sparse(int n) {
	/* This function takes as argument an integer n and returns a new
	2-dimensional array of doubles of shape n times n. */
	double ** A = new double *[n];
	for (int i = 0; i < n; i++) {
		A[i] = new double[n]();
	}
	return A;
}

void clear_memory(double ** & A, int n) {
	/* This function deletes 2 dimensional double arrays.*/
	for (int i = 0; i < n; i++) {
		delete[] A[i];
	}
	delete[] A;
}

double* zeros(int n) {
	/* This function takes as argument an integer n and returns a new 1-dimensional
	array of doubles of size n.*/
	double * v = new double[n]();
	return v;
}

class body {
	/*Represents a celestial body and its associated data.*/
public:
	double x = 0.0; double y = 0.0; double z = 0.0; //coordinates
	double vx = 0.0; double vy = 0.0; double vz = 0.0; //velocity
	double ax = 0.0; double ay = 0.0; double az = 0.0; //acceleration
	double m = 0.0; //mass
	std::string name;
};

using body_list = std::vector<body>;

class planetary_system {
	/* Represents a system of celestial bodies at a given time step.
	Contains functionality to update system to next time step, using either Euler forward or velocity Velvet.*/
private:
	double ** Fx; double ** Fy; double ** Fz; // forces in system
	double ** GMM; // mass product matrix Mi*Mj times G (gravitational constant)
	double * vxhalf; double * vyhalf; double * vzhalf; // half steps (in case of velocity velvet)

	// function for calculating (r cross v) angular momentum per mass of body i 
	double ang_mom_per_mass_squared(int i) {
		body body_i = bodies[i];
		double r_cross_v_x = (body_i.y*body_i.vz - body_i.z*body_i.vy);
		double r_cross_v_y = (body_i.z*body_i.vx - body_i.x*body_i.vz);
		double r_cross_v_z = (body_i.x*body_i.vy - body_i.y*body_i.vx);
		double l2 = r_cross_v_x * r_cross_v_x + r_cross_v_y * r_cross_v_y + r_cross_v_z * r_cross_v_z;
		return l2;
	}

	// function to calculate forces and find acceleration of celestial bodies at the given time step
	void update_acceleration() {
		// update forces
		for (int i = 0; i < n; i++) {
			body body_i = bodies[i];
			for (int j = 0; j < n; j++) {
				if (i != j && i > j) {
					body body_j = bodies[j];
					double xji = body_j.x - body_i.x;
					double zji = body_j.z - body_i.z;
					double yji = body_j.y - body_i.y;
					double r = sqrt(xji * xji + yji * yji + zji * zji);
					double r3 = pow(r, beta);
					Fx[i][j] = GMM[i][j] / r3 * xji;
					Fy[i][j] = GMM[i][j] / r3 * yji;
					Fz[i][j] = GMM[i][j] / r3 * zji;
					Fx[j][i] = -Fx[i][j];
					Fy[j][i] = -Fy[i][j];
					Fz[j][i] = -Fz[i][j];

					// adding optional relativistc correction to force
					if (relativistic == true) {
						double factor = three_by_c2 / (r*r);
						// body i
						Fx[i][j] *= 1 + ang_mom_per_mass_squared(i) * factor;
						Fy[i][j] *= 1 + ang_mom_per_mass_squared(i) * factor;
						Fz[i][j] *= 1 + ang_mom_per_mass_squared(i) * factor;

						// body j
						Fx[j][i] *= 1 + ang_mom_per_mass_squared(j) * factor;
						Fy[j][i] *= 1 + ang_mom_per_mass_squared(j) * factor;
						Fz[j][i] *= 1 + ang_mom_per_mass_squared(j) * factor;
					}
				}
			}
		}
		// calculate acceleration
		for (int i = 0; i < n; i++) {
			bodies[i].ax = 0.0; bodies[i].ay = 0.0; bodies[i].az = 0.0;
			for (int j = 0; j < n; j++) {
				double m = bodies[i].m;
				bodies[i].ax += Fx[i][j] / m;
				bodies[i].ay += Fy[i][j] / m;
				bodies[i].az += Fz[i][j] / m;
			}
		}
	}


public:
	// Physical properties of system
	int n; // number of bodies in system
	body_list bodies; // bodies in system (with associated position/velocity/mass data)

	// Physical quantities
	double time; // time
	double kinetic_energy; double potential_energy; // kinetic and potentetial energy
	double angular_momentum_x; double angular_momentum_y; double angular_momentum_z; // angular momentum in x, y, and z direction

	// Physical constants
	double beta = 3; // determines the power of 1/r^2. Default is 3 corresponding to 1/r^2
	double G = 1.985336e-29; // AU^3 kg^-1 years^-2 (gravitational constant)
	double c = 63197.791; // AU years^-1 (speed of light)
	double three_by_c2 = 3.0/(c*c);

	// Relativistic
	bool relativistic = false; // if set to true, then forces are found using relativistic term

	// constructor
	planetary_system(body_list data, double t) {
		time = t; bodies = data; n = data.size();
		Fx = sparse(n); Fy = sparse(n); Fz = sparse(n); GMM = sparse(n);
		vxhalf = zeros(n); vyhalf = zeros(n); vzhalf = zeros(n);

		//calculate mass products
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				GMM[i][j] = G*bodies[i].m*bodies[j].m;
			}
		}

		//initiate acceleration
		update_acceleration();
	}

	// destructor
	~planetary_system() {
		delete[] vxhalf; delete[] vyhalf; delete[] vzhalf; 
		clear_memory(Fx,n); clear_memory(Fy,n); clear_memory(Fz,n);
		clear_memory(GMM, n);
	}

	// function for calculating the next step after a time dt using Euler forward
	void next_step_eulerfd(double dt) {
		// find acceleration (given the position data)
		update_acceleration();

		// fix sun
		bodies[0].ax = 0; bodies[0].ay = 0; bodies[0].az = 0;
		bodies[0].vx = 0; bodies[0].vy = 0; bodies[0].vz = 0;

		for (int i = 0; i < n; i++) {
			// update velocity
			bodies[i].vx += bodies[i].ax*dt;
			bodies[i].vy += bodies[i].ay*dt;
			bodies[i].vz += bodies[i].az*dt;
			// update position
			bodies[i].x += bodies[i].vx*dt;
			bodies[i].y += bodies[i].vy*dt;
			bodies[i].z += bodies[i].vz*dt;
		}

		time += dt;
	}

	// function for calculating the next step after a time dt using velocity Verlet
	void next_step_velverlet(double dt) {
		for (int i = 0; i < n; i++) {
			// update vhalf
			vxhalf[i] = bodies[i].vx + 0.5*bodies[i].ax*dt;
			vyhalf[i] = bodies[i].vy + 0.5*bodies[i].ay*dt;
			vzhalf[i] = bodies[i].vz + 0.5*bodies[i].az*dt;

			// update position
			bodies[i].x += vxhalf[i] * dt;
			bodies[i].y += vyhalf[i] * dt;
			bodies[i].z += vzhalf[i] * dt;

		}
		// update acceleration (given the position data)
		update_acceleration();

		// fix sun
		bodies[0].ax = 0; bodies[0].ay = 0; bodies[0].az = 0;
		bodies[0].vx = 0; bodies[0].vy = 0; bodies[0].vz = 0;

		// update velocity
		for (int i = 0; i < n; i++) {
			bodies[i].vx = vxhalf[i] + 0.5*bodies[i].ax * dt;
			bodies[i].vy = vyhalf[i] + 0.5*bodies[i].ay * dt;
			bodies[i].vz = vzhalf[i] + 0.5*bodies[i].az * dt;
		}

		time += dt;
	}

	// function for calculating kinetic energy and storing the result as "kinetic_energy"
	void calculate_kinetic_energy() {
		double double_kinetic_energy_sum=0;
		for (int i = 0; i < n; i++) {
			body body_i = bodies[i];
			double v2 = body_i.vx*body_i.vx + body_i.vy*body_i.vy + body_i.vz*body_i.vz; //velocity squared
			double_kinetic_energy_sum += body_i.m*v2;
		}
		kinetic_energy = 0.5*double_kinetic_energy_sum;
	}

	// function for calculating the potential energy and storing the result as "potential_energy"
	void calculate_potential_energy() {
		double potential_energy_sum=0; 		
		for (int i = 0; i < n; i++) {
			body body_i = bodies[i];
			for (int j = 0; j < n; j++) {
				if (i != j && i>j) {
					body body_j = bodies[j];
					double xji = body_j.x - body_i.x;
					double zji = body_j.z - body_i.z;
					double yji = body_j.y - body_i.y;
					double r = sqrt(xji * xji + yji * yji + zji * zji);
					potential_energy_sum += -GMM[i][j]/r;
				}
			}
		}
		potential_energy = 2 * potential_energy_sum;
	}

	// function for calculating the total angular momentum
	void calculate_angular_momentum() {
		double r_cross_p_x = 0; double r_cross_p_y = 0; double r_cross_p_z = 0;
		for (int i=0; i < n; i++) {
			body body_i = bodies[i];
			r_cross_p_x += body_i.m*(body_i.y*body_i.vz - body_i.z*body_i.vy);
			r_cross_p_y += body_i.m*(body_i.z*body_i.vx - body_i.x*body_i.vz);
			r_cross_p_z += body_i.m*(body_i.x*body_i.vy - body_i.y*body_i.vx);
		}
		angular_momentum_x = r_cross_p_x;
		angular_momentum_y = r_cross_p_y;
		angular_momentum_z = r_cross_p_z;
	}
};

void initialize(body_list & bodies, std::string filename) {
	// function to initiate an instance of the class body_list, according to a
	// correctly formated file containing data from the solar system.
	std::ifstream inFile(filename);
	std::string line;

	while (std::getline(inFile, line))
	{
		std::stringstream templine(line);
		std::vector<std::string > words;
		std::string word;
		while (std::getline(templine, word, ' '))
		{
			words.push_back(word);
		}
		body celestical_body;
		celestical_body.name = words[0];
		celestical_body.m = std::stod(words[1]);
		celestical_body.x = std::stod(words[2]);
		celestical_body.y = std::stod(words[3]);
		celestical_body.z = std::stod(words[4]);
		celestical_body.vx = std::stod(words[5]) * 365;
		celestical_body.vy = std::stod(words[6]) * 365;
		celestical_body.vz = std::stod(words[7]) * 365;
		bodies.push_back(celestical_body);
	}
	inFile.close();
}

void write_conserved_quantities_to_file(std::ofstream &file, 
	body_list bodies, 
	std::string method_name, //either "euler" or "verlet"
	double number_years=1, 
	double dt=1.0/365) {
	/*Writes writes time (t), kinetic energy (K) and potential energy (U),
	angular momentum (Lx,Ly,Lz) to a CSV.*/

	file.precision(15);

	planetary_system bodies_system(bodies, 0.0);
	double t; // time
	double K; // kinetic energy
	double U; // potential energy
	double Lx; double Ly; double Lz; // angular momentum

	double day = 1.0/ 365;
	double ndays = 0.0;

	std::cout << "Beginning calculations for " << method_name << std::endl;
	while (bodies_system.time <= number_years) {
		// calculate energy and angular momentum
		if (bodies_system.time >= ndays) {
			bodies_system.calculate_kinetic_energy();
			bodies_system.calculate_potential_energy();
			bodies_system.calculate_angular_momentum();
			t = bodies_system.time;
			K = bodies_system.kinetic_energy;
			U = bodies_system.potential_energy;
			Lx = bodies_system.angular_momentum_x;
			Ly = bodies_system.angular_momentum_y;
			Lz = bodies_system.angular_momentum_z;
			file << t << ',' << K << ',' << U << ',' << Lx << ',' << Ly << ',' << Lz << std::endl;
			ndays += day;
		}
		if (method_name == "euler") {
			bodies_system.next_step_eulerfd(dt);
		}
		else if (method_name == "verlet") {
			bodies_system.next_step_velverlet(dt);
		}
	}
	std::cout << "Completed calculations for " << method_name << std::endl;
}

void unit_test_1(std::string method_name) {
	body_list earth_sun;
	initialize(earth_sun, "initial_conditions_earth_sun.txt");
	planetary_system earth_sun_system(earth_sun, 0.0);
	body earth_begin = earth_sun_system.bodies[1];
	double dt = 1.0 / (365 * 24);
	if (method_name == "euler") {
		while (earth_sun_system.time <= 1.0) {
			earth_sun_system.next_step_eulerfd(dt);
		};
	}
	else if (method_name == "verlet") {
		while (earth_sun_system.time <= 1.0) {
			earth_sun_system.next_step_velverlet(dt);
		};
	}
	body earth_final = earth_sun_system.bodies[1];
	if (abs(earth_begin.x - earth_final.x) < 1e-5 && abs(earth_begin.y - earth_final.y) < 1e-5) {
		std::cout << "Unit test " << method_name  <<" successful." << std::endl;
	}
}


void unit_test_converved_quantities(std::string method_name) {
	body_list earth_sun;
	initialize(earth_sun, "initial_conditions_earth_sun.txt");
	planetary_system earth_sun_system(earth_sun, 0.0);

	earth_sun_system.calculate_kinetic_energy();
	earth_sun_system.calculate_potential_energy();
	earth_sun_system.calculate_angular_momentum();

	double K_begin = earth_sun_system.kinetic_energy;
	double U_begin = earth_sun_system.potential_energy;
	double Lz_begin = earth_sun_system.angular_momentum_z;

	// loop
	double dt = 1.0 / (365 * 24);
	if (method_name == "euler") {
		while (earth_sun_system.time <= 1.0) {
			earth_sun_system.next_step_eulerfd(dt);
		};
	}
	else if (method_name == "verlet") {
		while (earth_sun_system.time <= 1.0) {
			earth_sun_system.next_step_velverlet(dt);
		};
	}

	earth_sun_system.calculate_kinetic_energy();
	earth_sun_system.calculate_potential_energy();
	earth_sun_system.calculate_angular_momentum();

	double K_end = earth_sun_system.kinetic_energy;
	double U_end = earth_sun_system.potential_energy;
	double Lz_end = earth_sun_system.angular_momentum_z;

	if (abs(K_end-K_begin)/abs(K_begin) < 1e-5 && abs(U_end - U_begin) / abs(U_begin) < 1e-5 && abs(Lz_end - Lz_begin) / abs(Lz_begin) < 1e-5) {
		std::cout << "Unit test " << method_name  << " successful." << std::endl;
	}
}

int main() {
	int N; 
	double number_years; double dt;
	double x; double y; double z;

	//======================================================================================================
	// Unit tests
	//======================================================================================================

	unit_test_1("euler");
	unit_test_1("verlet");
	unit_test_converved_quantities("euler");
	unit_test_converved_quantities("verlet");

	//======================================================================================================
	// Reading initial conditions from various files
	//======================================================================================================

	body_list earth_sun;
	initialize(earth_sun, "initial_conditions_earth_sun.txt");

	body_list earth_sun_jupiter;
	initialize(earth_sun_jupiter, "initial_conditions_earth_sun_jupiter.txt");

	body_list all_planets;
	initialize(all_planets, "initial_conditions_planets.txt");

	body_list all_planets_and_moons;
	initialize(all_planets_and_moons, "initial_conditions_planets_and_moons.txt");

	body_list mercury_sun;
	// read the file "initial_conditions_mercury_sun.txt" and insert its values into bodies
	initialize(mercury_sun, "initial_conditions_mercury_sun.txt");

	//=========================================================================================================
	//	Make x,y data for different values of dt of earth-sun system using verlet and euler
	//=========================================================================================================

	//std::vector<double> dt_values = { 1.0 / 12, 1.0 / 52, 1.0 / 365, 1.0 / (365 * 24), 1.0 / (365 * 24 * 4) ,1.0 / (365 * 24 * 60), 1.0 / (365 * 24 * 60 * 60), 1.0 / (365 * 24 * 60 * 60 * 10), (1.0 / (365 * 24 * 60 * 60 * 10))/10 }; // big wtf. 1.0 / (365 * 24 * 60 * 60 * 100) yields a negative number????
	//int length = size(dt_values);

	//std::cout << dt_values[7]  << std::endl;


	//// opening two files for each dt in dt_values, one for Euler and one for Verlet
	//std::vector<std::ofstream> data_euler(length);
	//std::vector<std::ofstream> data_verlet(length);
	//for (int i = 0; i < length; i++) {
	//	double dt = dt_values[i];
	//	data_euler[i].open("euler_dt_equals_" + std::to_string(i) + ".txt");
	//	data_verlet[i].open("verlet_dt_equals_" + std::to_string(i) + ".txt");
	//}

	////number_years = 100;
	//number_years = 1; 	

	//// writing euler data
	//std::cout << "starting Euler" << std::endl;

	//for (int i = 0; i < length; i++) {
	//	planetary_system earth_sun_system(earth_sun, 0.0);
	//	double dt = dt_values[i];
	//	std::cout << "started dt =" << dt << std::endl;
	//	double day = 1.0 / 365;
	//	double ndays = 0.0;
	//	while (earth_sun_system.time < number_years) {

	//		
	//		if (earth_sun_system.time >= ndays) { //print only data each day to file
	//			x = earth_sun_system.bodies[1].x;
	//			y = earth_sun_system.bodies[1].y;
	//			data_euler[i] << x << ',' << y << std::endl;
	//			ndays += day;
	//		}
	//		earth_sun_system.next_step_eulerfd(dt);
	//	}
	//	std::cout << "completed dt =" << dt << std::endl;
	//}

	//// writing verlet data
	//std::cout << "starting Verlet" << std::endl;
	//for (int i = 0; i < length; i++) {
	//	planetary_system earth_sun_system(earth_sun, 0.0);
	//	double dt = dt_values[i];
	//	std::cout << "started dt =" << dt << std::endl;
	//	double day = 1.0 / 365;
	//	double ndays = 0.0;
	//	while (earth_sun_system.time < number_years) { 
	//		

	//		if (earth_sun_system.time >= ndays) { //print only data each day to file
	//			x = earth_sun_system.bodies[1].x;
	//			y = earth_sun_system.bodies[1].y;
	//			data_verlet[i] << x << ',' << y << std::endl;
	//			ndays += day;
	//		}
	//		earth_sun_system.next_step_velverlet(dt);
	//	}
	//	std::cout << "completed dt =" << dt << std::endl;
	//}

	//// closing files
	//for (int i = 0; i < length; i++) {
	//	data_euler[i].close();
	//	data_verlet[i].close();
	//}

	//=========================================================================================================
	//	Time algorithm
	//=========================================================================================================

	//std::vector<body_list> test_body_lists = { earth_sun , earth_sun_jupiter, all_planets, all_planets_and_moons };
	//std::vector<std::string> plot_names = { "earth_sun", "earth_sun_jupiter", "all_planets", "all_planets_and_moons" };

	//std::vector<std::ofstream> times_euler(plot_names.size());
	//std::vector<std::ofstream> times_verlet(plot_names.size());

	//// opening files
	//std::cout <<"Opening files." << std::endl;
	//for (int i = 0; i < plot_names.size(); i++) {
	//	times_euler[i].open("euler_" + plot_names[i]+ "_time.txt");
	//	times_verlet[i].open("verlet_" + plot_names[i] + "_time.txt");
	//}

	//for (int i = 0; i < plot_names.size(); i++) {
	//	std::cout << "Calculating " << plot_names[i] << std::endl;
	//	using namespace std::chrono;
	//    
	//	body_list temp_body_list = test_body_lists[i];

	//	double timelimit = 1;
	//	double dt = 1.0 / (365*24*60);
	//	
	//	int number_laps = 10;
	//	for (int laps = 0; laps <= number_laps; laps++) {
	//		planetary_system temp_system_verlet(temp_body_list, 0.0);
	//		planetary_system temp_system_euler(temp_body_list, 0.0);

	//		// verlet
	//		high_resolution_clock::time_point t1 = high_resolution_clock::now();
	//		while (temp_system_verlet.time < timelimit) {
	//			temp_system_verlet.next_step_velverlet(dt);
	//		}
	//		high_resolution_clock::time_point t2 = high_resolution_clock::now();
	//		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	//		times_verlet[i] << time_span.count() << std::endl;

	//		// euler
	//		t1 = high_resolution_clock::now();
	//		while (temp_system_euler.time < timelimit) {
	//			temp_system_euler.next_step_eulerfd(dt);
	//		}
	//		t2 = high_resolution_clock::now();
	//		time_span = duration_cast<duration<double>>(t2 - t1);
	//		times_euler[i] << time_span.count() << std::endl;
	//	}

	//}

	//// closing files
	//std::cout << "Closing files." << std::endl;
	//for (int i = 0; i < plot_names.size(); i++) {
	//	times_euler[i].close();
	//	times_verlet[i].close();
	//}

	//=========================================================================================================
	//	Make energy,angular momentum data for set value of dt using verlet and euler
	//=========================================================================================================

	//std::ofstream E_L_data_euler("E_L_data_euler.txt");
	//std::ofstream E_L_data_verlet("E_L_data_verlet.txt");

	//write_conserved_quantities_to_file(E_L_data_euler,
	//	earth_sun,
	//	"euler",
	//	number_years=350, 
	//	dt = 1.0 / (365*24*60)); // 1 minute
	//write_conserved_quantities_to_file(E_L_data_verlet,
	//	earth_sun,
	//	"verlet",
	//	number_years = 350,
	//	dt = 1.0 / (365 * 24 * 60)); // 1 minute

	//E_L_data_euler.close();
	//E_L_data_verlet.close();

	//=========================================================================================================
	//	Make x,y data for all planets in the solar system
	//=========================================================================================================

	//planetary_system solar_system(earth_sun, 0.0);
	//planetary_system solar_system(all_planets, 0.0);
	//planetary_system solar_system(all_planets_and_moons, 0.0);
	//dt = 1.0 / (365*24*4);

	//// opening files for each body in the solar system, e.g., "verlet_earth.txt"
	//int n = solar_system.n;
	//std::vector<std::ofstream> myfiles(n);
	//for (int i = 0; i < n; i++) {
	//	myfiles[i].open("verlet_" + solar_system.bodies[i].name + ".txt");
	//	std::cout << "verlet_" + solar_system.bodies[i].name + ".txt" << std::endl;
	//}

	//// writing data to files
	//double day = 1.0 / 365;
	//double n_days = day;
	//while (solar_system.time <=3){
	//	if (solar_system.time> n_days){
	//		for (int i = 0; i < n; i++) {
	//			body currbody = solar_system.bodies[i];
	//			x = currbody.x;
	//			y = currbody.y;
	//			myfiles[i] << x << ',' << y << std::endl;
	//		}
	//		n_days += day;
	//	}
	//	solar_system.next_step_velverlet(dt);
	//}

	//// closing files
	//for (int i = 0; i < n; i++) {
	//	myfiles[i].close();
	//}

	//=========================================================================================================
	//	Study escape velocity/ different beta
	//=========================================================================================================

	//dt = 1.0 / (364 * 24); number_years = 10;
	//std::vector<double> s = { 0, 0.25, 0.5, 0.75, 1, 1.25}; std::vector<double> beta = { 3,3.25, 3.5, 3.75, 4 };
	//
	//
	//// writing data to files for different s and beta pairs
	//for (int i = 0; i < s.size(); i++) {
	//	body_list earth_sun;
	//	initialize(earth_sun, "initial_conditions_earth_sun.txt"); // load circular orbit
	//	earth_sun[1].vy = (1-s[i])*earth_sun[1].vy + s[i] * sqrt(2)*earth_sun[1].vy; // set new velocity
	//	for (int j = 0; j < beta.size(); j++) {
	//		std::string filename = "earth_escape_s_" + std::to_string(i) + "_beta_" + std::to_string(j) + ".txt";
	//		std::ofstream myfile(filename);
	//		planetary_system earth_sun_system(earth_sun, 0.0); // initiate system
	//		earth_sun_system.beta = beta[j];
	//		std::cout << "beginning plot with s=" << s[i] << " and beta=" << beta[j] - 1 << std::endl;
	//		while (earth_sun_system.time <= number_years) {
	//			x = earth_sun_system.bodies[1].x;
	//			y = earth_sun_system.bodies[1].y;
	//			myfile << x << ',' << y << std::endl;
	//			earth_sun_system.next_step_velverlet(dt);
	//		}
	//		myfile.close();
	//	}
	//}

	//=========================================================================================================
	//	 Study effect of Jupiters gravity on earth
	//=========================================================================================================

	std::vector<int> masses = { 1,10,100,1000 };
	
	for (int i = 0; i < size(masses); i++) {
		body_list earth_sun_jupiter;
		initialize(earth_sun_jupiter, "initial_conditions_earth_sun_jupiter.txt");
		earth_sun_jupiter[2].m = earth_sun_jupiter[2].m*masses[i];
		std::ofstream E_L_data_verlet_earth_jupiter("E_L_data_verlet_earth_jupiter_m_"+std::to_string(masses[i])+".txt");
		write_conserved_quantities_to_file(E_L_data_verlet_earth_jupiter,
			earth_sun_jupiter,
			"verlet",
			number_years = 350,
			dt = 1.0 / (365*24*60));
		E_L_data_verlet_earth_jupiter.close();
		dt = 1.0 / (365 * 24 * 60);
		planetary_system earth_sun_jupiter_system(earth_sun_jupiter, 0.0);
		std::ofstream x_y_earth_sun_jupiter_system_sun("x_y_earth_sun_jupiter_system_sun_m_" + std::to_string(masses[i]) + ".txt");
		std::ofstream x_y_earth_sun_jupiter_system_earth("x_y_earth_sun_jupiter_system_earth_m_" + std::to_string(masses[i]) + ".txt");
		std::ofstream x_y_earth_sun_jupiter_system_jupiter("x_y_earth_sun_jupiter_system_jupiter_m_" + std::to_string(masses[i]) + ".txt");
		double day = 1.0 / 365;
		double days = 0.0;
		while (earth_sun_jupiter_system.time <= 12) {
			if (earth_sun_jupiter_system.time>= days) {
				body sun= earth_sun_jupiter_system.bodies[0];
				body earth = earth_sun_jupiter_system.bodies[1];
				body jupiter = earth_sun_jupiter_system.bodies[2];
				x_y_earth_sun_jupiter_system_sun << sun.x << ',' << sun.y << std::endl;
				x_y_earth_sun_jupiter_system_earth << earth.x << ',' << earth.y << std::endl;
				x_y_earth_sun_jupiter_system_jupiter << jupiter.x << ',' << jupiter.y << std::endl;
				days += day;
			}
			earth_sun_jupiter_system.next_step_velverlet(dt);
		}
		x_y_earth_sun_jupiter_system_sun.close();
		x_y_earth_sun_jupiter_system_earth.close();
		x_y_earth_sun_jupiter_system_jupiter.close();
	}

	//=========================================================================================================
	//	 Study perihelion precession
	//=========================================================================================================

	
	//planetary_system mercury_sun_system(mercury_sun, 0.0);
	//mercury_sun_system.relativistic = true;
	//dt = 1.0 / (365*24*60*60); // 1 seconds

	//// opening initial conditions for mercury-sun system from file
	//std::ofstream prehelion_precession_file;
	//prehelion_precession_file.open("mercury_prehelion_precission.txt");
	//
	//// initialize previous step
	//double rprev = sqrt(mercury_sun_system.bodies[1].x*mercury_sun_system.bodies[1].x+ mercury_sun_system.bodies[1].y*mercury_sun_system.bodies[1].y);
	//mercury_sun_system.next_step_velverlet(dt);

	//// initialize current step
	//body mercury_curr = mercury_sun_system.bodies[1];
	//x = mercury_curr.x; y = mercury_curr.y;
	//double t = mercury_sun_system.time;
	//double rcurr = sqrt(mercury_curr.x*mercury_curr.x + mercury_curr.y*mercury_curr.y);

	//// writing data to files
	//while (mercury_sun_system.time <= 100.0){
	//	mercury_sun_system.next_step_velverlet(dt);
	//	body mercury = mercury_sun_system.bodies[1];
	//	double rnext = sqrt(mercury.x*mercury.x + mercury.y*mercury.y);
	//	if (rcurr < rnext && rcurr < rprev) {
	//		prehelion_precession_file << t << ','<< x << ',' << y << std::endl;
	//	}
	//	x = mercury.x; y = mercury.y;
	//	t = mercury_sun_system.time;
	//	rprev = rcurr;
	//	rcurr = rnext;
	//}

	//// closing files
	//prehelion_precession_file.close();

	//=========================================================================================================
	//	 success message if the problem runs
	//=========================================================================================================

	std::cout << "success!" << std::endl;
	std::cin.get();
}


