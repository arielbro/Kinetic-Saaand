using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kinetic_Saaand
{
    interface IParticlePairsSelector
    {
        IEnumerable<Tuple<Particle, Particle>> selectPairs(IEnumerable<Particle> particles);
    }
}
