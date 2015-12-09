//Written by Pirouz Nourian, researcher @ TU Delft, 3D Geoinformation in 2014. The work reported here has been funded by NLeSC Big Data Analytics in the Geo-Spatial Domain project (code: 027.013.703).
// the marching cube algorithm has been adapted from the c++ script by Paul Bourke and a VB.net code written by David Stasiuk 
using System;
using System.Drawing;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino;
using Rhino.Geometry;
using Grasshopper;

namespace TUD_RasterWorks_V_0_1
{
    public class Voxel
    {
        private double X;
        private double Y;
        private double Z;
        private double M;
        //= 1
        private double ResX;
        //= 1
        private double ResY;
        //= 1
        private double ResZ;
        //add a field for resultion (multi resolution octree and regular voxel grid)
        public Voxel(double XCo, double YCo, double ZCo, double Val)
        {
            this.X = XCo;
            this.Y = YCo;
            this.Z = ZCo;
            this.M = Val;
        }
        public Voxel(double XCo, double YCo, double ZCo, double Val, double XS, double YS, double ZS)
        {
            this.X = XCo;
            this.Y = YCo;
            this.Z = ZCo;
            this.M = Val;
            ResX = XS;
            ResY = YS;
            ResZ = ZS;
        }
        public Voxel(double XCo, double YCo, double ZCo, double Val, double Size)
        {
            this.X = XCo;
            this.Y = YCo;
            this.Z = ZCo;
            this.M = Val;
            ResX = Size;
            ResY = Size;
            ResZ = Size;
        }
        public Voxel(Point3d Vertex, double Value)
        {
            this.X = Vertex.X;
            this.Y = Vertex.Y;
            this.Z = Vertex.Z;
            this.M = Value;
        }
        public Voxel(Point3d Vertex, double Value, double XS, double YS, double ZS)
        {
            this.X = Vertex.X;
            this.Y = Vertex.Y;
            this.Z = Vertex.Z;
            this.M = Value;
            ResX = XS;
            ResY = YS;
            ResZ = ZS;
        }
        public Voxel(Point3d Vertex, double Value, double Size)
        {
            this.X = Vertex.X;
            this.Y = Vertex.Y;
            this.Z = Vertex.Z;
            this.M = Value;
            ResX = Size;
            ResY = Size;
            ResZ = Size;
        }
        public double Value
        {
            get { return this.M; }
            set { this.M = value; }
        }
        public Point3d vertex
        {
            get { return new Point3d(this.X, this.Y, this.Z); }
            set
            {
                this.X = value.X;
                this.Y = value.Y;
                this.Z = value.Z;
            }
        }
    }
    public class Raster3D
    {
        //Inherits Grasshopper.Kernel.Types.GH_Goo<>
        public Voxel[, ,] voxels;
        public List<Voxel> Voxellist;
        private Point3d[] points;
        private double[] Measures;
        private int XC;
        private int YC;
        private int ZC;
        //has to be corrected to make a voxel(,,) array as well. Otherwise will be of no use at all!
        public Raster3D(List<Point3d> Vertices, List<double> Values, int XCount, int YCount, int ZCount)
        {
            if (Vertices.Count != Values.Count)
            {
                //Print("the number of points and measures should be the same!!!")
                return;
            }
            int RefLength = XCount * YCount * ZCount;
            Voxel[] LocVoxels = new Voxel[RefLength];
            //As New List(Of Voxel)
            Voxel[, ,] TheVoxels = new Voxel[XCount, YCount, ZCount];
            int i = 0;
            int j = 0;
            int k = 0;
            for (int ind = 0; ind <= Vertices.Count - 1; ind++)
            {
                Voxel vx = new Voxel(Vertices[ind], Values[ind]);
                LocVoxels[ind] = (vx);
                i = ind % XCount;
                j = ind / XCount % YCount;
                k = ind / XCount * YCount % (XCount * YCount - 1);
                TheVoxels[i, j, k] = vx;
            }
            this.Voxellist = LocVoxels.ToList();
            this.points = Vertices.ToArray();
            this.Measures = Values.ToArray();
            this.XC = XCount;
            this.YC = YCount;
            this.ZC = ZCount;
            this.voxels = TheVoxels;
        }
        public Raster3D(Box Bounds, string FieldFunction, double Sx, double Sy, double Sz)
        {
            bool NoEval = false;
            if (FieldFunction == null)
            {
                NoEval = true;
            }
            Box BBox = Bounds;
            int i = 0;
            int j = 0;
            int k = 0;
            double u = 0;
            double v = 0;
            double w = 0;
            int m = 0;
            int n = 0;
            int l = 0;
            u = BBox.X.Length;
            v = BBox.Y.Length;
            w = BBox.Z.Length;
            m = Convert.ToInt32(Math.Floor(u / Sx));
            n = Convert.ToInt32(Math.Floor(v / Sy));
            l = Convert.ToInt32(Math.Floor(w / Sz));
            double spx = 0;
            double spy = 0;
            double spz = 0;
            //corrected sizes
            spx = u / m;
            spy = v / n;
            spz = w / l;
            Voxel[, ,] voxels = new Voxel[m + 1, n + 1, l + 1];
            int RefLength = (m + 1) * (n + 1) * (l + 1);
            Point3d[] vertices = new Point3d[RefLength];
            double[] values = new double[RefLength];
            string CorrectExpression = Grasshopper.Kernel.Expressions.GH_ExpressionSyntaxWriter.RewriteAll(FieldFunction);
            Grasshopper.Kernel.Expressions.GH_ExpressionParser Expr = new Grasshopper.Kernel.Expressions.GH_ExpressionParser();
            Expr.CacheSymbols(CorrectExpression);
            int Counter = 0;
            double Value = 0;
            Point3d VXL = default(Point3d);
            for (k = 0; k <= l; k++)
            {
                for (j = 0; j <= n; j++)
                {
                    for (i = 0; i <= m; i++)
                    {
                        VXL = BBox.PointAt(i * spx / u, j * spy / v, k * spz / w);
                        vertices[Counter] = (VXL);
                        if (NoEval)
                        {
                            Value = 0;
                        }
                        else
                        {
                            Expr.ClearVariables();
                            Expr.AddVariable("x", VXL.X);
                            Expr.AddVariable("y", VXL.Y);
                            Expr.AddVariable("z", VXL.Z);
                            Value = Expr.Evaluate(Expr.CachedSymbols())._Double;
                        }
                        values[Counter] = (Value);
                        voxels[i, j, k] = new Voxel(VXL, Value);
                        Counter = Counter + 1;
                    }
                }
            }
            this.XC = m + 1;
            this.YC = n + 1;
            this.ZC = l + 1;
            this.points = vertices;
            this.Measures = values;
        }
        public Voxel[, ,] VolumetricPixels
        {
            get { return this.voxels; }
            set { this.voxels = value; }
        }
        public Point3d[] VoxelPoints
        {
            get { return this.points; }
            set { this.points = value; }
        }
        public double[] VoxelMeasures
        {
            get { return this.Measures; }
            set { this.Measures = value; }
        }
        public int XCount
        {
            get { return this.XC; }
        }
        public int YCount
        {
            get { return this.YC; }
        }
        public int ZCount
        {
            get { return this.ZC; }
        }
        public List<Voxel> VoxList
        {
            get { return this.Voxellist; }
        }

        public Int32[] intEdgeTable = 
      {
      0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,      0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
      0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,      0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
      0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,      0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
      0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,      0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
      0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,      0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
      0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,      0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
      0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,      0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
      0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,      0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
      0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,      0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
      0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,      0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
      0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,      0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
      0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,      0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
      0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,      0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
      0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,      0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
      0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,      0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
      0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,      0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
      };

        Int32[,] triTable =
      {
      {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
      {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
      {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
      {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
      {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
      {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
      {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
      {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
      {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
      {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
      {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
      {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
      {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
      {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
      {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
      {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
      {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
      {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
      {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
      {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
      {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
      {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
      {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
      {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
      {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
      {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
      {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
      {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
      {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
      {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
      {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
      {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
      {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
      {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
      {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
      {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
      {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
      {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
      {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
      {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
      {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
      {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
      {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
      {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
      {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
      {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
      {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
      {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
      {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
      {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
      {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
      {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
      {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
      {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
      {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
      {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
      {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
      {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
      {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
      {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
      {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
      {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
      {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
      {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
      {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
      {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
      {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
      {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
      {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
      {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
      {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
      {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
      {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
      {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
      {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
      {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
      {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
      {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
      {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
      {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
      {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
      {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
      {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
      {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
      {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
      {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
      {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
      };

        public Mesh Field_to_IsoSurface(double ts)
        {
            Raster3D R3D = this;

            double[][] GMs = null;
            // can be converted to a 2D array
            Point3d[][] GBs = Grid_to_Boxes(Plane.WorldXY, R3D, ref GMs);
            Mesh IsoSurface = new Mesh();
            for (int GI = 0; GI <= GBs.Length - 1; GI++)
            {
                Mesh facet = local_cube(GBs[GI], GMs[GI], ts);
                IsoSurface.Append(facet);
            }
            IsoSurface.Vertices.CombineIdentical(true, true);
            //    IsoSurface.Normals.ComputeNormals()
            //    IsoSurface.UnifyNormals()
            return IsoSurface;
        }
        public Point3d[][] Grid_to_Boxes(Plane BasePlane, Raster3D R3D, ref double[][] BoxMs)
        {
            Voxel[, ,] Voxels = R3D.voxels;
            int XCount = 0;
            int YCount = 0;
            int ZCount = 0;
            XCount = R3D.XCount;
            YCount = R3D.YCount;
            ZCount = R3D.ZCount;
            int BCount = (XCount - 1) * (YCount - 1) * (ZCount - 1);
            Point3d[][] GBs = new Point3d[BCount][];
            double[][] GMs = new double[BCount][];
            int i = 0;
            int j = 0;
            int k = 0;
            int Counter = 0;
            for (k = 0; k <= ZCount - 2; k++)
            {
                for (j = 0; j <= YCount - 2; j++)
                {
                    for (i = 0; i <= XCount - 2; i++)
                    {
                        Point3d[] points = new Point3d[8];
                        // = New Point3d(7){}
                        double[] values = new double[8];
                        // = New Double(7){}
                        points[0] = R3D.voxels[i, j, k].vertex;
                        values[0] = R3D.voxels[i, j, k].Value;
                        points[1] = R3D.voxels[i + 1, j, k].vertex;
                        values[1] = R3D.voxels[i + 1, j, k].Value;
                        points[2] = R3D.voxels[i + 1, j + 1, k].vertex;
                        values[2] = R3D.voxels[i + 1, j + 1, k].Value;
                        points[3] = R3D.voxels[i, j + 1, k].vertex;
                        values[3] = R3D.voxels[i, j + 1, k].Value;
                        points[4] = R3D.voxels[i, j, k + 1].vertex;
                        values[4] = R3D.voxels[i, j, k + 1].Value;
                        points[5] = R3D.voxels[i + 1, j, k + 1].vertex;
                        values[5] = R3D.voxels[i + 1, j, k + 1].Value;
                        points[6] = R3D.voxels[i + 1, j + 1, k + 1].vertex;
                        values[6] = R3D.voxels[i + 1, j + 1, k + 1].Value;
                        points[7] = R3D.voxels[i, j + 1, k + 1].vertex;
                        values[7] = R3D.voxels[i, j + 1, k + 1].Value;
                        GBs[Counter] = points;
                        GMs[Counter] = values;
                        Counter += 1;
                    }
                }
            }
            BoxMs = GMs;
            return GBs;
        }
        private Mesh local_cube(Point3d[] c, double[] c_v, double ts)
        {
            //marching cube algorithm
            //adapted from the c++ script by Paul Bourke and a VB.NET code written by David Stasiuk 
            //http://paulbourke.net/geometry/polygonise/, http://www.grasshopper3d.com/profiles/blogs/marching-cubes-curve-wrapping-more-metaballs :
    

            Int32 cubeindex = default(Int32);
            if (c_v[0] <= ts)
                cubeindex = cubeindex | 1;
            if (c_v[1] <= ts)
                cubeindex = cubeindex | 2;
            if (c_v[2] <= ts)
                cubeindex = cubeindex | 4;
            if (c_v[3] <= ts)
                cubeindex = cubeindex | 8;
            if (c_v[4] <= ts)
                cubeindex = cubeindex | 16;
            if (c_v[5] <= ts)
                cubeindex = cubeindex | 32;
            if (c_v[6] <= ts)
                cubeindex = cubeindex | 64;
            if (c_v[7] <= ts)
                cubeindex = cubeindex | 128;

            Point3d[] vertlist = new Point3d[12];
            // edge middle points, there are 12 edges in every cube

            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 1))
                vertlist[0] = interp_vertex(c[0], c[1], c_v[0], c_v[1], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 2))
                vertlist[1] = interp_vertex(c[1], c[2], c_v[1], c_v[2], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 4))
                vertlist[2] = interp_vertex(c[2], c[3], c_v[2], c_v[3], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 8))
                vertlist[3] = interp_vertex(c[3], c[0], c_v[3], c_v[0], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 16))
                vertlist[4] = interp_vertex(c[4], c[5], c_v[4], c_v[5], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 32))
                vertlist[5] = interp_vertex(c[5], c[6], c_v[5], c_v[6], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 64))
                vertlist[6] = interp_vertex(c[6], c[7], c_v[6], c_v[7], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 128))
                vertlist[7] = interp_vertex(c[7], c[4], c_v[7], c_v[4], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 256))
                vertlist[8] = interp_vertex(c[0], c[4], c_v[0], c_v[4], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 512))
                vertlist[9] = interp_vertex(c[1], c[5], c_v[1], c_v[5], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 1024))
                vertlist[10] = interp_vertex(c[2], c[6], c_v[2], c_v[6], ts);
            if (Convert.ToBoolean(intEdgeTable[cubeindex] & 2048))
                vertlist[11] = interp_vertex(c[3], c[7], c_v[3], c_v[7], ts);

            Mesh box_mesh = new Mesh();

            for (Int32 i = 0; i <= 16; i += 3)
            {
                if (triTable[cubeindex, i] == -1)
                    break; // TODO: might not be correct. Was : Exit For
                box_mesh.Vertices.Add(vertlist[triTable[cubeindex, i]]);
                box_mesh.Vertices.Add(vertlist[triTable[cubeindex, i + 1]]);
                box_mesh.Vertices.Add(vertlist[triTable[cubeindex, i + 2]]);
                box_mesh.Faces.AddFace(new MeshFace(i + 2, i + 1, i));
            }

            return box_mesh;

        }

        private Point3d interp_vertex(Point3d p1, Point3d p2, double v1, double v2, double ts)
        {
            return new Point3d(p1 + (ts - v1) / (v2 - v1) * (p2 - p1));
        }
        public bool refine_mesh(ref Mesh msh)
        {
            msh.Vertices.CombineIdentical(true, true);
            msh.FaceNormals.ComputeFaceNormals();
            msh.Normals.ComputeNormals();
            msh.UnifyNormals();
            return true;
        }

        public Mesh Raster3D_to_Raster2D(List<System.Drawing.Color> Colors, double param, int xyz)
        {
            //Dim R3D As OTB_3DGIS.Raster3d = Nothing
            Raster3D R3D = this;
            //CType(R3D_O, OTB_3DGIS.Raster3D)
            Voxel[, ,] Slice = R3D.voxels;
            int i = 0;
            int j = 0;
            int k = 0;
            int m = 0;
            int n = 0;
            int l = 0;
            m = R3D.XCount - 1;
            n = R3D.YCount - 1;
            l = R3D.ZCount - 1;
            Point3d[,] PR2D_XP = new Point3d[n + 1, l + 1];
            Point3d[,] PR2D_YP = new Point3d[l + 1, m + 1];
            Point3d[,] PR2D_ZP = new Point3d[m + 1, n + 1];
            List<Color> Cols = new List<Color>();
            List<Point3d> PixelPoints = new List<Point3d>();
            Mesh ZM = new Mesh();
            if (xyz == 0)
            {
                Interval PDomain = new Interval(0, l);
                int KP = Convert.ToInt32(Math.Round(PDomain.ParameterAt(param)));
                int counter = KP * (n + 1) * (l + 1);
                for (k = 0; k <= l; k++)
                {
                    for (j = 0; j <= n; j++)
                    {
                        Voxel ThatVoxel = R3D.voxels[KP, j, k];
                        Point3d PixelPoint = ThatVoxel.vertex;
                        PR2D_XP[j, k] = PixelPoint;
                        PixelPoints.Add(PixelPoint);
                        Cols.Add(Colors[counter]);
                        counter = counter + 1;
                    }
                }
                ZM = MeshFromPoints(PixelPoints.ToList(), Cols, n + 1, l + 1);
            }
            else if (xyz == 1)
            {
                Interval PDomain = new Interval(0, l);
                int KP = Convert.ToInt32(Math.Round(PDomain.ParameterAt(param)));
                int counter = KP * (l + 1) * (m + 1);
                for (k = 0; k <= l; k++)
                {
                    for (i = 0; i <= m; i++)
                    {
                        Voxel ThatVoxel = R3D.voxels[i, KP, k];
                        Point3d PixelPoint = ThatVoxel.vertex;
                        PR2D_YP[k, i] = PixelPoint;
                        PixelPoints.Add(PixelPoint);
                        Cols.Add(Colors[counter]);
                        counter = counter + 1;
                    }
                }
                ZM = MeshFromPoints(PixelPoints.ToList(), Cols, l + 1, m + 1);
            }
            else if (xyz == 2)
            {
                Interval PDomain = new Interval(0, l);
                int KP = Convert.ToInt32(Math.Round(PDomain.ParameterAt(param)));
                int counter = KP * (m + 1) * (n + 1);
                for (j = 0; j <= n; j++)
                {
                    for (i = 0; i <= m; i++)
                    {
                        Voxel ThatVoxel = R3D.voxels[i, j, KP];
                        Point3d PixelPoint = ThatVoxel.vertex;
                        PR2D_ZP[i, j] = PixelPoint;
                        PixelPoints.Add(PixelPoint);
                        Cols.Add(Colors[counter]);
                        counter = counter + 1;
                    }
                }
                ZM = MeshFromPoints(PixelPoints.ToList(), Cols, m + 1, n + 1);
            }
            else
            {
                return null;
            }
            return ZM;
        }
        public Mesh MeshFromPoints(List<Point3d> P, List<Color> C, int U, int V)
        {
            //, C As list(Of Color)
            Mesh M = new Mesh();
            M.Vertices.AddVertices(P);
            if ((C != null))
                M.VertexColors.AppendColors(C.ToArray());
            for (int j = 0; j <= V - 2; j++)
            {
                for (int i = 0; i <= U - 2; i++)
                {
                    int n0 = 0;
                    int n1 = 0;
                    int n2 = 0;
                    int n3 = 0;
                    n0 = j * U + i;
                    n1 = n0 + 1;
                    n2 = n0 + U;
                    n3 = n1 + U;
                    MeshFace face = new MeshFace(n0, n1, n3, n2);
                    M.Faces.AddFace(face);
                }
            }
            return M;
        }
    }
}
