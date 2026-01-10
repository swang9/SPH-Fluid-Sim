#define GL_SILENCE_DEPRECATION
#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <vector>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

// solver parameters
const static Vector2d G(0.f, -10.f);   // external force (2×1 vector of type double)
const static float REST_DENS = 300.f;  // rest density
const static float GAS_CONST = 2000.f; // const for equation of state
const static float H = 16.f;		   // kernel radius
const static float HSQ = H * H;		   // radius^2 for optimization
const static float MASS = 2.5f;		   // assume all particles have the same mass
const static float VISC = 200.f;	   // viscosity constant
const static float DT = 0.0007f;	   // integration timestep

// SPH smoothing kernels
const static float POLY6 = 4.f / (M_PI * pow(H, 8.f)); // poly6 
const static float SPIKY_GRAD = -10.f / (M_PI * pow(H, 5.f)); // spiky
const static float VISC_LAP = 40.f / (M_PI * pow(H, 5.f)); // viscosity

// simulation parameters
const static float EPS = H; // boundary epsilon
const static float BOUND_DAMPING = -0.5f;
const static int DAM_PARTICLES = 500;


// rendering projection parameters
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static double VIEW_WIDTH = 1.5 * 800.f;
const static double VIEW_HEIGHT = 1.5 * 600.f;



// PARTICLE STRUCT
struct Particle
{
    // position (x), velocity (v), force (f), density (rho), pressure (p)
    // uses Eigen library - Vector2d: 2×1 vector of type double
	Particle(float _x, float _y) : x(_x, _y), v(0.f, 0.f), f(0.f, 0.f), rho(0), p(0.f) {}
	Vector2d x, v, f;
	float rho, p;
};

// solver data
static vector<Particle> particles;

// RENDERING SYSTEMS
void InitGL(void){
    // initializes GL: sets background color, point smoothing, 
    glClearColor(0.9f, 0.9f, 0.9f, 1);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(H / 2.f);
	glMatrixMode(GL_PROJECTION);
}
void Render(void){
    // clears GL and sets up coord system, particle color
    // loops through particles array and places vertex at x/y
    glClear(GL_COLOR_BUFFER_BIT);

	glLoadIdentity();
	glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);

	glColor4f(0.2f, 0.6f, 1.f, 1);
	glBegin(GL_POINTS);
	for (auto &p : particles)
	{
		glVertex2f(p.x(0), p.x(1));
	}
	glEnd();

	glutSwapBuffers();
}

// SOLVER
void InitSPH(void){
    cout << "initializing dam break with " << DAM_PARTICLES << " particles" << endl;
	for (float y = EPS; y < VIEW_HEIGHT - EPS * 2.f; y += H)
	{
		for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H)
		{
			if (particles.size() < DAM_PARTICLES)
			{
                // add randomness
				float jitter = static_cast<float>(arc4random()) / static_cast<float>(RAND_MAX);
				particles.push_back(Particle(x + jitter, y));
			}
			else
			{
				return;
			}
		}
	}
}
void ComputeDensityPressure(void){
    // loops through all particles, summing their density contributions
    for (auto &pi : particles)
	{
		pi.rho = 0.f;
		for (auto &pj : particles)
		{
			Vector2d rij = pj.x - pi.x;
			float r2 = rij.squaredNorm();

			if (r2 < HSQ)
			{
				pi.rho += MASS * POLY6 * pow(HSQ - r2, 3.f);
			}
		}
        // ideal gas law
		pi.p = GAS_CONST * (pi.rho - REST_DENS);
	}
}
void ComputeForces(void){
    // loops through all particles, summing their pressure and viscosity
    for (auto &pi : particles)
	{
		Vector2d fpress(0.f, 0.f);
		Vector2d fvisc(0.f, 0.f);
		for (auto &pj : particles)
		{
			if (&pi == &pj)
			{
				continue;
			}

			Vector2d rij = pj.x - pi.x;
			float r = rij.norm();

			if (r < H)
			{
				// compute pressure force contribution
				fpress += -rij.normalized() * MASS * (pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD * pow(H - r, 3.f);
				// compute viscosity force contribution
				fvisc += VISC * MASS * (pj.v - pi.v) / pj.rho * VISC_LAP * (H - r);
			}
		}
		Vector2d fgrav = G * MASS / pi.rho;
		pi.f = fpress + fvisc + fgrav;
	}
}
void Integrate(void){
    // loops through all particles and applies Euler's method
    for (auto &p : particles)
    {
        // forward Euler integration
        p.v += DT * p.f / p.rho;
        p.x += DT * p.v;

        // enforce boundary conditions
        if (p.x(0) - EPS < 0.f)
        {
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = EPS;
        }
        if (p.x(0) + EPS > VIEW_WIDTH)
        {
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = VIEW_WIDTH - EPS;
        }
        if (p.x(1) - EPS < 0.f)
        {
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = EPS;
        }
        if (p.x(1) + EPS > VIEW_HEIGHT)
        {
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = VIEW_HEIGHT - EPS;
        }
    }
}
void Update(void){
    // called between frame renders to step sim state forward
    // according to SPH, first compute densities, then compute pressure
    // then use densities and pressure to calculate the forces acting on each particle, which are used by Integrate()
    ComputeDensityPressure();
	ComputeForces();
	Integrate();

	glutPostRedisplay();
}

// MAIN
int main(int argc, char **argv)
{
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInit(&argc, argv);
	glutCreateWindow("Müller SPH");

	// sets render to generate new frames as often as possible
	glutDisplayFunc(Render);
	// sets update to move sim forward
	glutIdleFunc(Update);

	// initilizes rendering environement with particles
    InitGL();
	InitSPH();

	glutMainLoop();
	return 0;
}