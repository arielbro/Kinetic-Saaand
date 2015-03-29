using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kinetic_Saaand
{
    public class Particle : IComparable<Particle>
    {
        private double mass;
        private int dim;
        private double[] position;
        private double[] velocity;
        private double[] totalForce;

        public Particle(double mass, double[] startPosition, double[] startVelocity=null)
        {
            this.mass = mass;
            this.dim = startPosition.Length;
            this.position = startPosition;
            this.velocity = (double[])startVelocity.Clone() ?? new double[dim];
        }

        /// <summary>
        /// Returns a vector of size dim holding the particle's position.
        /// </summary>
        public double[] getPosition()
        {
            return (double[])position.Clone();
        }


        /// <summary>
        /// Returns a vector of dim holding the particle's velocity.
        /// </summary>
        public double[] getVelocity()
        {
            return (double[])velocity.Clone();
        }

        /// <summary>
        /// Sets the position of the particle.
        /// </summary>
        /// <param name="vector">vector of size dim holding the new particle position</param>
        public void setPosition(double[] newPosition)
        {
            position = (double[])newPosition.Clone();
        }

        /// <summary>
        /// Sets the velocity of the particle.
        /// </summary>
        /// <param name="vector">vector of size dim holding the new particle velocity</param>
        public void setVelocity(double[] newVelocity)
        {
            velocity = (double[])newVelocity.Clone();
        }


        public double getMass()
        {
            return mass;
        }

        public int CompareTo(Particle other)
        {
            if (other == null)
                return 1;
            //implement order lexicographically - first coordinate has most weight, then other, etc.
            double[] thisPosition = this.getPosition();
            double[] otherPosition = other.getPosition();
            for (int coordinate = 0; coordinate < dim; coordinate++)
            {
                if (thisPosition[coordinate] != otherPosition[coordinate])
                    if(thisPosition[coordinate] > otherPosition[coordinate])
                        return 1;
                    else
                        return -1;
            }
            //particles with same position are declared equal, for that matter
            return 0;
        }

        public static bool operator <(Particle e1, Particle e2)
        {
            return e1.CompareTo(e2) < 0;
        }

        public static bool operator >(Particle e1, Particle e2)
        {
            return e1.CompareTo(e2) > 0;
        }
    }
}


