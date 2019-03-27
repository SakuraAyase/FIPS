#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#define pi 3.1415926535
#define E  2.71828182845904523536

using namespace std;

double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}
class Particle
{
public:
	//double fitness;
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;
	vector<double>Pm;
	//double pBestFitness;
	Particle() {}

	Particle(vector<double> position, vector<double>velocity, vector<double>best_position, double best_fitness)
	{
		this->position = position;
		this->velocity = velocity;
		this->pBest = best_position;
		//this->pBestFitness = best_fitness;
	}
};

bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

class PSO
{
public:

	PSO(int dim, int m, int Tmax, double max, double min, double c, 
		 double dt, double percent,double X)
	{
		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->c = c;
		this->dt = dt;
		this->percent = percent;
		this->X = X;
		particles.resize(m);

	}


	double fitnessFunction(vector<double> pos)
	{
		double result = 0.0;
		for (int i = 0; i < dim; i++)
		{
			result += pow(pos[i], 2);
		}

		/*
		double temp1 = 0.0;
		double temp2 = 0.0;
		for (int i = 0; i < dim; i++)
		{
			double x = particle.position[i];
			temp1 += pow(x, 2);
			temp2 += cos(2 * pi*x);
		}

		result = -20 * exp(-0.2*sqrt(temp1 / dim)) + 20 + E - exp(temp2 / dim);
		*/
		return result;
	}

	void initialParticles(int i)
	{
		particles[i].position.resize(dim);
		particles[i].velocity.resize(dim);
		particles[i].pBest.resize(dim);
		particles[i].Pm.resize(dim);
		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);
			particles[i].position[j] = randDouble(this->min, this->max);
			particles[i].velocity[j] = randDouble(-range, range);
			particles[i].pBest[j] = particles[i].position[j];
		}
		//particles[i].fitness = fitnessFunction(particles[i].position);
		//particles[i].pBestFitness = fitnessFunction(particles[i].pBest);
	}

	void initialAllParticles()
	{

		for (int i = 0; i < m; i++)
		{
			initialParticles(i);
		}
	}

	void inertiaWeight()
	{
		//w = randDouble(0.4, 0.6);
		double t = T / ((double)Tmax);
		w = 0.9;
	}

	void updateParticle(int i)
	{
		int left = i - 1;
		if (left < 0)
			left = m - 1;
		int right = (i + 1) % m;
		double wleft = fitnessFunction(particles[left].pBest);
		double wright = fitnessFunction(particles[right].pBest);
		double wself = fitnessFunction(particles[i].pBest);
		double k1 = randDouble(0, c / 2);
		double k2 = randDouble(0, c / 2);
		double k3 = randDouble(0, c / 2);
		for (int j = 0; j < dim; j++) 
		{
			particles[i].Pm[j] = (particles[i].pBest[j] * k3 * wself + 
				particles[left].pBest[j] * k1 * wleft +
				particles[right].pBest[j] * k2*wright) / (wright*k2 + wleft * k1 + k3 * wself);
		}

		//cout << particles[i].Pm[0] << " " << particles[i].position[0]<< endl;
		
		for (int j = 0; j < dim; j++)
		{
			double last_position = particles[i].position[j];
			double range = percent * (max - min);

			particles[i].velocity[j] = X*(particles[i].velocity[j] +
				c * (particles[i].Pm[j] - particles[i].position[j]));
			particles[i].position[j] += dt * particles[i].velocity[j];

			if (particles[i].velocity[j] > range)
				particles[i].velocity[j] = range;
			if (particles[i].velocity[j] < -range)
				particles[i].velocity[j] = -range;

			if (particles[i].position[j] > max)
				particles[i].position[j] = max;
			if (particles[i].position[j] < min)
				particles[i].position[j] = min;

		}
		//particles[i].fitness = fitnessFunction(particles[i].position);

		if (fitnessFunction(particles[i].position) < fitnessFunction(particles[i].pBest))
		{
			//particles[i].pBestFitness = particles[i].fitness;

			particles[i].pBest = particles[i].position;
		}

	}


	void updateAllParticles()
	{
		//inertiaWeight();
		for (int i = 0; i < m; i++)
		{
			updateParticle(i);
		}
		T++;
	}

	double getFitness()
	{
		int index = 0;
		for (int i = 0; i < m; i++)
		{
			//cout << fitnessFunction(particles[i].position) << endl;
			if (fitnessFunction(particles[i].pBest) < fitnessFunction(particles[index].pBest))
				index = i;
		}
		//cout << index << endl;
		return fitnessFunction(particles[index].pBest);
	}
private:
	int dim;
	int m;//number of instances

	int T;
	int Tmax;

	double w;
	double max;
	double min;
	double c;

	double X;

	double dt;//时间步长
	double percent;


	vector<Particle> particles;


};

void run(vector<double>& result1)
{
	int dim = 30;
	int m = 20;
	int Tmax = 20000;
	double max = 100;
	double min = -100;
	double c = 4.1;
	double dt = 1.0;
	double percent = 0.2;
	double X = 0.729;

	PSO pso = PSO(dim, m, Tmax, max, min, c, dt, percent,X);
	pso.initialAllParticles();

	vector<double>fitness;
	fitness.push_back(pso.getFitness());

	for (int i = 0; i < Tmax; i++)
	{
		pso.updateAllParticles();
		cout << ":";
		//fitness.push_back(pso.getFitness());
		fitness.push_back(pso.getFitness());
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << pso.getFitness() << endl;
	}

	result1 = fitness;
}

int main()
{

	int times = 5;
	int interval = 10;
	vector<double> result1;

	run(result1);

	for (int i = 1; i < times; i++)
	{
		vector<double> result1_temp;
		run(result1_temp);
		for (int j = 0; j < result1_temp.size(); j++)
		{
			result1[j] += result1_temp[j];
		}
	}
	for (int j = 0; j < result1.size(); j++)
	{
		result1[j] /= times;
	}

	for (int j = 0; j < result1.size(); j++)
	{
		if (j%interval == 0)
			cout << result1[j] << " ";
	}

	system("pause");
}