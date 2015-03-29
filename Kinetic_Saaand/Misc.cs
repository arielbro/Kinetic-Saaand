using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kinetic_Saaand
{
    static class Misc
    {
        public static IEnumerable<double[]> CartesianProduct(params double[][] inputs)
        {
            return (IEnumerable<double[]>)inputs.Aggregate(
                (IEnumerable<object[]>)new object[][] { new object[0] },
                (soFar, input) =>
                    from prevProductItem in soFar
                    from item in input
                    select prevProductItem.Concat(new object[] { item }).ToArray());
        }

    }

    public class DoubleArrayEqualityComparer : IEqualityComparer<double[]>
    {
        public bool Equals(double[] x, double[] y)
        {
            if (x.Length != y.Length)
            {
                return false;
            }
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] != y[i])
                {
                    return false;
                }
            }
            return true;
        }

        public int GetHashCode(double[] obj)
        {
            int result = 17;
            for (int i = 0; i < obj.Length; i++)
            {
                unchecked
                {
                    result = result * 23 + obj[i].GetHashCode();
                }
            }
            return result;
        }
    }

}
