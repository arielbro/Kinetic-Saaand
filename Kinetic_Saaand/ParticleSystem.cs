using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using WpfApplication1.Properties;

namespace Kinetic_Saaand
{
    //based on http://www.cs.cmu.edu/~baraff/pbm/particles.pdf
    class ParticleSystem
    {
        private List<Particle> particles; //Should figure out how to properly access particles as readonly (perhaps protected methods on them)
        private int dim;
        private Tuple<double, double>[] boundingBox;
        public event EventHandler TimeTick;
        
        public ParticleSystem(int dim, Tuple<double, double>[] boundingBox)
        {
            this.dim = dim;
            this.particles = new List<Particle>();
            this.boundingBox = boundingBox;
            
        }

        public void createParticle(double mass, double[] position, double[] velocity=null)
        {
            Debug.Assert(position.Length == dim);
            Particle particle = new Particle(mass, position, velocity);
            particles.Add(particle);
        }

        private double[] calculateGravity(Particle particle, int axis=1)
        {
            //constant pull from a hypothetical Earth
            double[] force = new double[dim];
            force[axis] = -Settings.Default.gravityConstant * particle.getMass();
            return force;
        }

        private double[] calculateDrag(Particle particle, double deltaT, double kDrag=0.1)
        {
            //calculate viscous drag - the effect of media (air/fluid) resistance to a moving particle. 
            //formula (by wolfram alpha) - 1/2 * kDrag * mass_density * area * speed ^ 2
            double[] velocity = particle.getVelocity();
            double[] force = new double[dim];
            for (int axis = 0; axis < dim; axis++)
            {
                double calculatedDragForce = -0.5 * Math.Sign(velocity[axis]) * kDrag * velocity[axis]*velocity[axis];
                //double maxDragAForce = -velocity[axis] * particle.getMass() / deltaT;
                force[axis] = calculatedDragForce;
            }
            return force;
        }

        private double[] calculateBoundingBoxForce(Particle particle, double[] actingForces, double deltaT, 
                                                   double frictionCoef=0.1)
        {
            /// <summary>
            /// Calculate the foces the bounding box exerts on a particle.
            /// for each dimension, there are two bounding hyperplanes. For each one, check if the particle is at contact/beyond it.
            /// If the particle has any velocity in the hyperplane's direction, negate it, with the coefficient of restitution
            /// defining how inelastic the collision is.  
            /// If the particle does not have velocity in this direction, it is in contact.
            /// Negate the component of the force acting on the particle in the hyperplane's direction, and enact linear drag force in the
            /// tangent part of the velocity. The force is proportional to the pressure on the surface (normal part of the particle force),
            /// the friction coefficient, and the particles tangent velocity, but will not be stronger than the force needed to halt the 
            /// particle.
            /// </summary>
            double[] boundingBoxForces = new double[dim];
            double[] particlePosition = particle.getPosition();
            double[] particleVelocity = particle.getVelocity();
            double particleMass = particle.getMass();
            for(int axis=0; axis<dim; axis++)
            {
                for(int side=-1; side<2; side=side+2)//in each axis there are two bounding hyperplanes
                {
                    if((side==-1 && particlePosition[axis] <= boundingBox[axis].Item1) || 
                        (side==1 && particlePosition[axis] >= boundingBox[axis].Item2))
                    {
                        if(side*particleVelocity[axis] >= Settings.Default.epsilon) //collision
                        {
                            boundingBoxForces[axis] -= particleVelocity[axis] * (1 + Settings.Default.kRestitute) * particleMass /
                                                       deltaT;
                            Console.WriteLine("Collision! axis: " + axis + " side: " + side + " axis velocity: " + particleVelocity[axis]
                                 + " axis position: " + particlePosition[axis]);
                        }
                        else if (Math.Abs(particleVelocity[axis]) < Settings.Default.epsilon) //contact
                        {
                            Console.WriteLine("Contact! axis: " + axis + " side: " + side + " axis velocity: " + particleVelocity[axis]
                                 + " axis position: " + particlePosition[axis]);                            //negate residual speed in that axis
                            boundingBoxForces[axis] += side * particleVelocity[axis] * particleMass /
                                                       deltaT;
                            //negate force normal to the hyperplane (pointing outwards)
                            if (side * actingForces[axis] > 0)
                                boundingBoxForces[axis] -= actingForces[axis];
                            //enact linear friction drag force (acts on each tangent velocity component)
                            for(int tangentAxis=0; tangentAxis<dim; tangentAxis++)
                            {
                                if(tangentAxis == axis)
                                    continue;
                                double calculatedDragForce = frictionCoef * -particleVelocity[tangentAxis] * 
                                    Math.Abs(boundingBoxForces[axis]);
                                double maxDragAForce = -particleVelocity[tangentAxis] * particleMass / deltaT;
                                boundingBoxForces[tangentAxis] += Math.Abs(calculatedDragForce) > Math.Abs(maxDragAForce) ?
                                                                  maxDragAForce : calculatedDragForce;
                            }
                        }
                        else //due to 
                        {
                            Console.WriteLine("Escape? axis: " + axis + " side: " + side + " axis velocity: " + particleVelocity[axis]
                                 + " axis position: " + particlePosition[axis]);
                        }
                    }
                }
            }
            return boundingBoxForces;
        }

        private double[][] calculateForces(double deltaT)
        {
            double[][] forces = new double[particles.Count][];
            for(int particleIndex=0; particleIndex<particles.Count; particleIndex++)
            {
                double[] forceAccumulator = new double[dim];
                double[] tempForce = calculateGravity(particles[particleIndex]);
                for (int i = 0; i < dim; i++)
                    forceAccumulator[i] += tempForce[i];
                tempForce = calculateDrag(particles[particleIndex], deltaT);
                for (int i = 0; i < dim; i++)
                    forceAccumulator[i] += tempForce[i];
                tempForce = calculateBoundingBoxForce(particles[particleIndex], forceAccumulator, deltaT);
                for (int i = 0; i < dim; i++)
                    forceAccumulator[i] += tempForce[i];
                forces[particleIndex] = forceAccumulator;
            }
            return forces;
        }

        /// <summary>
        /// Calculates and returns a matrix with dimensions (number of particles X 2*dim), with each row
        /// representing a particle's position and velocity.
        /// </summary>
        private double[][] getParticleStates()
        {
            //Should we just keep it all in one big scripty matrix for performance, and dump the attampt at OO?
            double[][] states = new double[particles.Count][];
            for(int i=0; i < particles.Count; i++)
            {
                Particle particle = particles[i];
                states[i] = particle.getPosition().Concat(particle.getVelocity()).ToArray();
            }
            return states;
        }

        /// <summary>
        /// Sets the position and velocity of all particles in the system.
        /// </summary>
        /// <param name="vector">matrix of size (number of particles X 2*dim), 
        /// with each row representing a particle's position and velocity
        /// </param>
        private void setParticleStates(double[][] states)
        {
            Debug.Assert(states.Length == particles.Count);
            for (int i = 0; i<states.Length ; i++)
            {
                Particle particle = particles[i];
                particle.setPosition(states[i].Take(dim).ToArray());
                particle.setVelocity(states[i].Skip(dim).ToArray());
            }
        }

        /// <summary>
        /// Uses the Euler method to approximate the states of the particles, given their current
        /// states and the forces acting on them.
        /// </summary>
        /// <param name="currentStates">matrix of dimensions (number of particles X 2*dim), representing the
        /// position and velocity of particles in the system</param>
        /// <param name="forces">matrix of dimensions (number of particles X 2*dim), representing the forces
        /// currently acting on theparticles in the system</param>
        /// <param name="deltaT">the length of the time interval after which new particle states are approximated
        /// </param>
        private double[][] performEulerStep(double[][] currentStates, double[][] forces, double deltaT)
        {
            Debug.Assert(currentStates.Length == particles.Count && forces.Length == particles.Count);
            Console.WriteLine("Tick");
            double[][] newStates = new double[particles.Count][];
            for(int i=0; i < particles.Count; i++)
            {
                Debug.Assert(currentStates[i].Length == 2*dim && forces[i].Length == dim);
                newStates[i] = new double[2*dim];
                //update positions
                for(int j = 0; j < dim; j++)
                    newStates[i][j] = currentStates[i][j] + currentStates[i][j + dim] * deltaT;
                //update velocities
                for (int j = 0; j < dim; j++)
                    newStates[i][j + dim] = currentStates[i][j + dim] + forces[i][j] * deltaT / particles[i].getMass();
            }
            return newStates;
        }

        /// <summary>
        /// Updates the state of the system by calculating current forces on particles, using a differential 
        /// solver to approximate the new states after applying forces and considering velocities.
        /// </summary>
        /// <param name="deltaT">the length of the time interval after which new particle states are approximated
        /// </param>
        public Task tick_time(double deltaT)
        {
            double[][] currentStates = getParticleStates();
            double[][] forces = calculateForces(deltaT);
            double[][] newStates = performEulerStep(currentStates, forces, deltaT);
            setParticleStates(newStates);
            return null;
        }
        public List<Particle> getParticles()
        {
            return particles.Select(item => item).ToList();//cloning
        }

        public double interpolateCollisionDeltaT(double particleAxisPosition, double particleAxisVelocity, double obstacleAxisPosition)
        {
            if(Math.Sign(obstacleAxisPosition - particleAxisPosition) != Math.Sign(particleAxisVelocity))
                return Double.NaN; //no collision, particle headed the other way (or at rest)
            return (obstacleAxisPosition - particleAxisPosition)/particleAxisVelocity;
        }
    }
}
