#include <iostream>
#include <vector> 
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

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
	/* This function deletes 2 dimensional double arrays.
	*/
	for (int i = 0; i < n; i++) {
		delete[] A[i];
	}
	delete[] A;
}

double* zeros(int n) {
	/* This function takes as argument an integer n and returns a new 1-dimensional
	array of doubles of size n.
	*/
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
	double ** MM; // mass product matrix Mi*Mj
	double * vxhalf; double * vyhalf; double * vzhalf; // half steps (in case of velocity velvet)

	// function to calculate forces and find acceleration of celestial bodies at the given time step
	void update_acceleration() {
		body body_i; body body_j;
		double xji; double yji; double zji;
		double r; double r3; double m;
		// update forces
		for (int i = 0; i < n; i++) {
			body_i = bodies[i];
			for (int j = 0; j < n; j++) {
				if (i != j && i > j) {
					body_j = bodies[j];
					xji = body_j.x - body_i.x;
					zji = body_j.z - body_i.z;
					yji = body_j.y - body_i.y;
					r = sqrt(xji * xji + yji * yji + zji * zji);
					r3 = pow(r, beta);
					Fx[i][j] = G * MM[i][j] / r3 * xji;
					Fy[i][j] = G * MM[i][j] / r3 * yji;
					Fz[i][j] = G * MM[i][j] / r3 * zji;
					Fx[j][i] = -Fx[i][j];
					Fy[j][i] = -Fy[i][j];
					Fz[j][i] = -Fz[i][j];
				}
			}
		}
		// calculate acceleration
		for (int i = 0; i < n; i++) {
			bodies[i].ax = 0.0; bodies[i].ay = 0.0; bodies[i].az = 0.0;
			for (int j = 0; j < n; j++) {
				m = bodies[i].m;
				bodies[i].ax += Fx[i][j] / m;
				bodies[i].ay += Fy[i][j] / m;
				bodies[i].az += Fz[i][j] / m;
			}
		}
	}


public:
	double time; // time
	body_list bodies; // bodies in system
	int n; // number of bodies in system
	double kinetic_energy; double potential_energy; // kinetic and potentetial energy
	double angular_momentum_x; double angular_momentum_y; double angular_momentum_z; // angular momentum in x, y, and z direction

	double beta = 3;
	double G = 1.985336e-29; // AU^3 kg^-1 years^-2 (gravitational constant)

	// constructor
	planetary_system(body_list data, double t) {
		time = t; bodies = data; n = data.size();
		Fx = sparse(n); Fy = sparse(n); Fz = sparse(n); MM = sparse(n);
		vxhalf = zeros(n); vyhalf = zeros(n); vzhalf = zeros(n);

		//calculate mass products
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				MM[i][j] = bodies[i].m*bodies[j].m;
			}
		}

		//initiate acceleration
		update_acceleration();
	}

	// destructor
	~planetary_system() {
		delete[] vxhalf; delete[] vyhalf; delete[] vzhalf; 
		clear_memory(Fx,n); clear_memory(Fy,n); clear_memory(Fz,n);
		clear_memory(MM, n);
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
		double half_kinetic_energy_sum=0; 
		double v2; //velocity squared
		for (int i = 0; i < n; i++) {
			v2 = bodies[i].vx*bodies[i].vx + bodies[i].vy*bodies[i].vy + bodies[i].vz*bodies[i].vz;
			half_kinetic_energy_sum += bodies[i].m*v2;
		}
		kinetic_energy = 0.5*half_kinetic_energy_sum;
	}

	// function for calculating the potential energy and storing the result as "potential_energy"
	void calculate_potential_energy() {
		body body_i; body body_j;
		double xji; double yji; double zji; double r2;
		double potential_energy_sum=0; 		
		for (int i = 0; i < n; i++) {
			body_i = bodies[i];
			for (int j = 0; j < n; j++) {
				if (i != j && i>j) {
					body_j = bodies[j];
					xji = body_j.x - body_i.x;
					zji = body_j.z - body_i.z;
					yji = body_j.y - body_i.y;
					r2 = xji * xji + yji * yji + zji * zji;
					potential_energy_sum += -MM[i][j]/r2;
				}
			}
		}
		potential_energy = 2 * G * potential_energy_sum;
	}

	// function for calculating the angular momentum
	void calculate_angular_momentum() {
		double r_cross_p_x = 0; double r_cross_p_y = 0; double r_cross_p_z = 0;
		body body_i;
		for (int i=0; i < n; i++) {
			body_i = bodies[i];
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
	std::string method_name,
	double number_years=1, 
	double dt=1.0/365) {

	planetary_system bodies_system(bodies, 0.0);
	double t; // time
	double K; // kinetic energy
	double U; // potential energy
	double Lx; double Ly; double Lz; // angular momentum
	while (bodies_system.time <= number_years) {
		// calculate energy and angular momentum
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
		if (method_name == "euler") {
			bodies_system.next_step_eulerfd(dt);
		}
		else if (method_name == "verlet") {
			bodies_system.next_step_velverlet(dt);
		}
	}
}

int main() {
	int N; 
	double number_years; double dt;
	double x; double y; double z;

	body_list earth_sun;
	// read the file "initial_conditions_earth_sun.txt" and insert its values into bodies
	initialize(earth_sun, "initial_conditions_earth_sun.txt");

	body_list all_planets_and_moons;
	// read the file "initial_conditions.txt" and insert its values into bodies
	initialize(all_planets_and_moons, "initial_conditions.txt");

	//=========================================================================================================
	//	Make x,y data for different values of dt of earth-sun system using verlet and euler
	//=========================================================================================================

	//std::vector<double> dt_values = { 1.0 / 12, 1.0 / 52, 1.0 / 365, 1.0 / (365 * 24) };//, 1.0 / (365 * 24 * 4), 1.0 / (365 * 24 * 60)};

	//// opening two files for each dt in dt_values, one for Euler and one for Verlet
	//int length = size(dt_values);
	//std::vector<std::ofstream> data_euler(length);
	//std::vector<std::ofstream> data_verlet(length);
	//for (int i = 0; i < length; i++) {
	//	double dt = dt_values[i];
	//	data_euler[i].open("euler_dt_equals_" + std::to_string(i) + ".txt");
	//	data_verlet[i].open("verlet_dt_equals_" + std::to_string(i) + ".txt");
	//}


	//number_years = 100;

	//// writing euler data
	//std::cout << "starting Euler" << std::endl;
	//for (int i = 0; i < length; i++) {
	//	planetary_system earth_sun_system(earth_sun, 0.0);
	//	dt = dt_values[i];
	//	std::cout << "started dt =" << dt << std::endl;
	//	while (earth_sun_system.time <= number_years) {
	//		x = earth_sun_system.bodies[1].x;
	//		y = earth_sun_system.bodies[1].y;
	//		data_euler[i] << x << ',' << y << std::endl;
	//		earth_sun_system.next_step_eulerfd(dt);
	//	}
	//	std::cout << "completed dt =" << dt << std::endl;
	//}

	//// writing verlet data
	//std::cout << "starting Verlet" << std::endl;
	//for (int i = 0; i < length; i++) {
	//	planetary_system earth_sun_system(earth_sun, 0.0);
	//	dt = dt_values[i];
	//	std::cout << "started dt =" << dt << std::endl;
	//	while (earth_sun_system.time <= number_years) {
	//		x = earth_sun_system.bodies[1].x;
	//		y = earth_sun_system.bodies[1].y;
	//		data_verlet[i] << x << ',' << y << std::endl;
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
	//	Make energy,angular momentum data for set value of dt using verlet and euler
	//=========================================================================================================

	//std::ofstream E_L_data_euler("E_L_data_euler.txt");
	//std::ofstream E_L_data_verlet("E_L_data_verlet.txt");

	//write_conserved_quantities_to_file(E_L_data_euler,
	//	earth_sun,
	//	"euler",
	//	number_years=350,
	//	dt = 1.0 / 365);
	//write_conserved_quantities_to_file(E_L_data_verlet,
	//	earth_sun,
	//	"verlet",
	//	number_years = 350,
	//	dt = 1.0 / 365);

	//E_L_data_euler.close();
	//E_L_data_verlet.close();

	//=========================================================================================================
	//	Make x,y data for all planet in the solar system
	//=========================================================================================================

	//planetary_system solar_system(earth_sun, 0.0);
	//solar_system.beta = 3.0;
	//dt = 1.0 / (365*24);

	//// opening files for each body in the solar system, e.g., "verlet_earth.txt"
	//int n = solar_system.n;
	//std::vector<std::ofstream> myfiles(n);
	//for (int i = 0; i < n; i++) {
	//	myfiles[i].open("verlet_" + solar_system.bodies[i].name + ".txt");
	//	std::cout << "verlet_" + solar_system.bodies[i].name + ".txt" << std::endl;
	//}

	//// writing data to files
	//while (solar_system.time <= 3.0){
	//	for (int i = 0; i < n; i++) {
	//		x = solar_system.bodies[i].x;
	//		y = solar_system.bodies[i].y;
	//		myfiles[i] << x << ',' << y << std::endl;
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
	//	body_list earth_sun_escape;
	//	initialize(earth_sun_escape, "initial_conditions_earth_sun.txt"); // load circular orbit
	//	earth_sun_escape[1].vy = (1-s[i])*earth_sun_escape[1].vy + s[i] * sqrt(2)*earth_sun_escape[1].vy; // set new velocity
	//	for (int j = 0; j < beta.size(); j++) {
	//		std::string filename = "earth_escape_s_" + std::to_string(i) + "_beta_" + std::to_string(j) + ".txt";
	//		std::ofstream myfile(filename);
	//		planetary_system earth_sun_system(earth_sun_escape, 0.0); // initiate system
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
			dt = 1.0 / 365);
		E_L_data_verlet_earth_jupiter.close();
		dt = 1.0 / 365;
		planetary_system earth_sun_jupiter_system(earth_sun_jupiter, 0.0);
		std::ofstream x_y_earth_sun_jupiter_system_earth("x_y_earth_sun_jupiter_system_earth_m_" + std::to_string(masses[i]) + ".txt");
		std::ofstream x_y_earth_sun_jupiter_system_jupiter("x_y_earth_sun_jupiter_system_jupiter_m_" + std::to_string(masses[i]) + ".txt");
		while (earth_sun_jupiter_system.time <= 12) {
			body earth = earth_sun_jupiter_system.bodies[1];
			body jupiter = earth_sun_jupiter_system.bodies[2];
			x_y_earth_sun_jupiter_system_earth << earth.x << ',' << earth.y << std::endl;
			x_y_earth_sun_jupiter_system_jupiter << jupiter.x << ',' << jupiter.y << std::endl;
			earth_sun_jupiter_system.next_step_velverlet(dt);
		}
		x_y_earth_sun_jupiter_system_earth.close();
		x_y_earth_sun_jupiter_system_jupiter.close();
	}

	//=========================================================================================================
	//	 success message if the problem runs
	//=========================================================================================================

	std::cout << "success!" << std::endl;
	std::cin.get();
}






//class body_list {
//	/*Represents a list of bodies and their associated data*/
//			//int n; //number of bodies in list
//	body* bodies;
//	int n;
//public:
//	// constructor
//	body_list(int size) {
//		n = size;
//		bodies = new body[n];
//	}
//	// return number of bodies
//	int size() {
//		return n;
//	}
//	// clear memory
//	void clear() {
//		delete[]bodies;
//	}
//	// access i'th body in bodies
//	body& operator[](int i) {
//		return bodies[i];
//	}
//};
















