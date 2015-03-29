using System;
using System.Collections.Generic;
using System.Configuration;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;
using WpfApplication1.Properties;

namespace Kinetic_Saaand
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private ParticleSystem particleSystem;
        Dictionary<Particle, GeometryModel3D> particleToGraphicMap;
        volatile bool isParticleSystemBusy = false;

        public MainWindow()
        {
            InitializeComponent();
            Tuple<double, double>[] boundingBox = new Tuple<double, double>[3];
            boundingBox[0] = new Tuple<double, double>(-35, 35);
            boundingBox[1] = new Tuple<double, double>(-10, 30);
            boundingBox[2] = new Tuple<double, double>(-30, 30);
            particleSystem = new ParticleSystem(3, boundingBox);
            particleToGraphicMap = new Dictionary<Particle, GeometryModel3D>();

            //create bounding box sides
            for(int axis=0; axis<3; axis++)
                for (int side = -1; side < 2; side++)
                {
                    continue;
                    //MeshGeometry3D mesh = new MeshGeometry3D();
                    //double axisValue = side > 0 ? boundingBox[axis].Item2 : boundingBox[axis].Item1;
                    //take all combinations of min/max of other axes with the axis value.
                    //mesh.Positions.Add();
                }

            //create particles
            Random random = new Random();
            for (int i = 0; i < 1; i++)
            {
                double[] startPosition = new double[] 
                {
                    random.NextDouble() * 40 - 20, 
                    random.NextDouble()*10 + 20,
                    random.NextDouble() * 40 - 20 
                };
                double[] startVelocity = new double[]
                {
                    random.NextDouble() * 30 - 15,
                    random.NextDouble() * 20 - 20,
                    random.NextDouble() * 30 - 15
                };
                double mass = random.NextDouble()*5 + 0.4;
                particleSystem.createParticle(mass, startPosition, startVelocity);
            }

            //create corresponding 3D models - note that all models start at origin.
            foreach (Particle particle in particleSystem.getParticles())
            {
                GeometryModel3D sphereGeometryModel = new GeometryModel3D();
                byte[] rgb = new byte[3];
                random.NextBytes(rgb);
                sphereGeometryModel.Material = new DiffuseMaterial(new SolidColorBrush(Color.FromRgb(rgb[0],rgb[1],rgb[2])));
                MeshGeometry3D geometry = new MeshGeometry3D();
                Sphere3D sphere = new Sphere3D(random.NextDouble()*0 + 0.15, Settings.Default.particleResolution);
                geometry.Positions = sphere.points;
                geometry.TriangleIndices = sphere.triangleIndices;
                sphereGeometryModel.Geometry = geometry;
                particleToGraphicMap.Add(particle, sphereGeometryModel);
            }

            //load models to viewport
            foreach (GeometryModel3D model in particleToGraphicMap.Values)
                model3DGroup.Children.Add(model);

            //add some fourth dimension magic (i.e. start the clock tick)
            System.Windows.Threading.DispatcherTimer dispatcherTimer = new System.Windows.Threading.DispatcherTimer();
            dispatcherTimer.Tick += new EventHandler(dispatcherTimer_Tick);
            dispatcherTimer.Interval = new TimeSpan(0, 0, 0, 0, Settings.Default.viewportRefreshRateMilliSec);
            dispatcherTimer.Start();
        }

        private async void dispatcherTimer_Tick(object sender, EventArgs e)
        {
            Task task = await Task.Factory.StartNew( () => particleSystem.tick_time(Settings.Default.deltaT));
            updatePositions();
        }

        private void updatePositions()
        {
            foreach (Particle particle in particleToGraphicMap.Keys)
            {
                double[] position = particle.getPosition();
                particleToGraphicMap[particle].Transform = new TranslateTransform3D(position[0], position[1], position[2]);
            }
        }
    }
}
