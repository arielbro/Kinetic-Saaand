using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kinetic_Saaand
{

    class DistanceParticlePairsSelector : IParticlePairsSelector
    {
        private int dim;
        private double maxDistance;

        public DistanceParticlePairsSelector(int dim, double maxDistance)
        {
            this.dim = dim;
            this.maxDistance = maxDistance;
        }


        //
        public IEnumerable<Tuple<Particle, Particle>> selectPairs(IEnumerable<Particle> particles)
        {
            List<Particle> particlesList = particles.ToList(); //for multiple iterations
            Tuple<double[], double[]> bounds = getCoordinateBounds(particlesList);
            double[] minBounds = bounds.Item1;
            double[] maxBounds = bounds.Item2;
            
            //create a grid with boxes of side length maxDistance, and divide the points to bins representing
            //the grid boxes. That guarentees that points in non-touching bins are of distance of atleast maxDistance.
            //By only generating non-empty bins, the number of bins is bounded by the number of points. The time complexity
            //of selecting all pairs of points from same/neighboring grids is O(M^2 * d * 2k), when M is the maximum number
            //of points in a grid, d is the number of non-empty grids, and k is the dimension. This is bounded by O(n^2).
            //If the points are somewhat evenly distributed, this makes for a big reduction in running time.
            
            //prepare bins
            Dictionary<double[], List<Particle>> gridsByPosition = new Dictionary<double[], List<Particle>>();
            foreach(Particle particle in particles)
            {
                double[] binLowerCorner = new double[dim]; //given a particle's position, the corresponding bin needs to be identified.
                double[] particlePosition = particle.getPosition();
                for(int i=0; i<dim; i++)
                {
                    double boxIndex = Math.Floor((particlePosition[i] - minBounds[i]) / maxDistance);
                    binLowerCorner[i] = minBounds[i] + maxDistance * boxIndex;
                }
                if (!gridsByPosition.ContainsKey(binLowerCorner))
                    gridsByPosition.Add(binLowerCorner, new List<Particle>());
                gridsByPosition[binLowerCorner].Add(particle);
            }
            //select pairs from bins
            HashSet<Tuple<Particle, Particle>> pairs = new HashSet<Tuple<Particle,Particle>>();
            HashSet<double[]> alreadyVisited = new HashSet<double[]>();
            foreach (double[] position in gridsByPosition.Keys)
            {
                double[][] adjacentCoordinates = new double[dim][];
                for(int coordinate =0 ; coordinate < dim; coordinate ++)
                    adjacentCoordinates[coordinate] = new double[]{
                        position[coordinate] - maxDistance, position[coordinate], position[coordinate] + maxDistance};
                foreach(double[] neighborPosition in Misc.CartesianProduct(adjacentCoordinates))
                {
                    if(alreadyVisited.Contains(neighborPosition))
                        continue;
                    if(neighborPosition == position)
                    {
                                    //pairs in the same bin
                    foreach(Particle particle1 in gridsByPosition[position])
                        foreach(Particle particle2 in gridsByPosition[neighborPosition])
                        if(particle2 < particle1) //don't repeat pairs (x,y == y,x)
                            pairs.Add(new Tuple<Particle, Particle>(particle1, particle2));
                    }
                    else //pairs in adjacent bins
                        foreach(Particle particle1 in gridsByPosition[position])
                            foreach(Particle particle2 in gridsByPosition[neighborPosition])
                                pairs.Add(new Tuple<Particle, Particle>(particle1, particle2));
          
                }
                alreadyVisited.Add(position);
            }
            return pairs;
        }

        private Tuple<double[], double[]> getCoordinateBounds(List<Particle> particles)
        {
            double[] minBounds = new double[dim];
            double[] maxBounds = new double[dim];
            minBounds = maxBounds = particles[0].getPosition();
            foreach(Particle particle in particles)
            {
                double[] position = particle.getPosition();
                for(int i=0; i<dim; i++)
                {
                    if (position[i] < minBounds[i])
                        minBounds[i] = position[i];
                    if (position[i] > maxBounds[i])
                        maxBounds[i] = position[i];
                }
            }
            return new Tuple<double[], double[]>(minBounds, maxBounds);
        }


    }
}
