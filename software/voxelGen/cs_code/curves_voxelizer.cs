//Written by Pirouz Nourian, researcher @ TU Delft, 3D Geoinformation in 2014. 
//The work reported here has been funded by NLeSC Big Data Analytics in the Geo-Spatial
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino;
using Rhino.Geometry;

namespace TUD_TopoVoxelizer_Curves_CSharp
{
    public class Class1
    {
        public List<Point3d> VoxelizePolyline(int CO, PolylineCurve Crv, double Sx, double Sy, double Sz)
        {
            //A = MeshToTriangles(M)
            Vector3d VSize = new Vector3d(Sx, Sy, Sz);
            BoundingBox RBBox = Crv.GetBoundingBox(false);
            BoundingBox ZBBox = RBBox_to_ZBBox(RBBox, VSize);
            //A = ZBBox
            List<Point3d> gridpoints = BBoxToVoxels(ZBBox, VSize);
            //B = Gridpoints
            List<Mesh> Intargets = new List<Mesh>();
            List<Point3d> Voxels = new List<Point3d>();
            int[] Facets = null;
            Mesh Intarget = null;
            List<Point3d> Crv_NearPoints = NearPoints(gridpoints, Crv, VSize.Length / 2);
            //print(Crv_NearPoints.count)
            foreach (Point3d gridpoint in Crv_NearPoints)
            {
                if (CO == 26)
                {
                    Intarget = LineTarget26_Connected(gridpoint, VSize);
                }
                else if (CO == 6)
                {
                    Intarget = LineTarget6_Connected(gridpoint, VSize);
                }
                else
                {
                    throw new Exception("connectivity target undefined!");
                }
                if (Rhino.Geometry.Intersect.Intersection.MeshPolyline(Intarget, Crv,out Facets).Length > 0)
                {
                    Voxels.Add(gridpoint);
                }
            }
            //B = Intargets
            return Voxels;
        }
        public Mesh[] MeshToTriangles(Mesh M)
        {
            M.Faces.ConvertQuadsToTriangles();
            Mesh[] Meshes = null;
            M.Unweld(0, false);
            Meshes = M.ExplodeAtUnweldedEdges();
            return Meshes;
        }
        public BoundingBox RBBox_to_ZBBox(BoundingBox RBbox, Vector3d VSize)
        {
            Point3d ZMinP = MinBoundRP_to_ZP(RBbox.Min, VSize);
            Point3d ZMaxP = MaxBoundRP_to_ZP(RBbox.Max, VSize);
            BoundingBox ZBBox = new BoundingBox(ZMinP, ZMaxP);
            return ZBBox;
        }
        public Point3d RPtoZP(Point3d RPoint, Vector3d VSize)
        {
            Point3d ZPOint = new Point3d((Math.Floor(RPoint.X / VSize.X) + 0.5) * VSize.X, (Math.Floor(RPoint.Y / VSize.Y) + 0.5) * VSize.Y, (Math.Floor(RPoint.Z) / VSize.Z + 0.5) * VSize.Z);
            return ZPOint;
        }
        private Point3d MaxBoundRP_to_ZP(Point3d RPoint, Vector3d VSize)
        {
            double x = 0;            double y = 0;            double z = 0;
            double u = 0;            double v = 0;            double w = 0;
            
            double ZP_x = 0;            double ZP_y = 0;            double ZP_z = 0;
            
            x = RPoint.X;            y = RPoint.Y;            z = RPoint.Z;
            
            u = VSize.X;            v = VSize.Y;            w = VSize.Z;
            
            ZP_x = Math.Ceiling(x / u) * u;            ZP_y = Math.Ceiling(y / v) * v;            ZP_z = Math.Ceiling(z / w) * w;
            Point3d ZPoint = new Point3d(ZP_x, ZP_y, ZP_z);
            return ZPoint;
        }
        private Point3d MinBoundRP_to_ZP(Point3d RPoint, Vector3d VSize)
        {
            double x = 0;            double y = 0;            double z = 0;           
            double u = 0;            double v = 0;            double w = 0;
            double ZP_x = 0;            double ZP_y = 0;            double ZP_z = 0;
            x = RPoint.X;            y = RPoint.Y;            z = RPoint.Z;
            u = VSize.X;            v = VSize.Y;            w = VSize.Z;
            ZP_x = Math.Floor(x / u) * u;            ZP_y = Math.Floor(y / v) * v;            ZP_z = Math.Floor(z / w) * w;
            Point3d ZPoint = new Point3d(ZP_x, ZP_y, ZP_z);
            return ZPoint;
        }
        public List<Point3d> BBoxToVoxels(BoundingBox ZBBox, Vector3d Vsize)
        {
            int i = 0;            int j = 0;            int k = 0;
            int Imax = 0;            int Jmax = 0;            int Kmax = 0;
            Imax = Convert.ToInt32(Math.Abs(ZBBox.Diagonal.X / Vsize.X) - 1);
            Jmax = Convert.ToInt32(Math.Abs(ZBBox.Diagonal.Y / Vsize.Y) - 1);
            Kmax = Convert.ToInt32(Math.Abs(ZBBox.Diagonal.Z / Vsize.Z) - 1);
            //Dim VoxelPoints(IMax,JMax,KMax) As Point3d
            List<Point3d> VoxelPoints = new List<Point3d>();
            for (k = 0; k <= Kmax; k++)
            {
                for (j = 0; j <= Jmax; j++)
                {
                    for (i = 0; i <= Imax; i++)
                    {
                        Point3d RelPoint = new Point3d(i * Vsize.X, j * Vsize.Y, k * Vsize.Z);
                        //VoxelPoints(i, j, k) = (RelPoint + New Vector3d(ZBbox.Min) + New Vector3d(VSize(0), VSize(1), VSize(2)) * 0.5)
                        VoxelPoints.Add((RelPoint + new Vector3d(ZBBox.Min) + new Vector3d(Vsize[0], Vsize[1], Vsize[2]) * 0.5));
                    }
                }
            }
            return VoxelPoints;
        }
        public Mesh LineTarget26_Connected(Point3d VoxelPoint, Vector3d Vsize)
        {
            Mesh IntersectionTarget = new Mesh();
            Point3d[] Vertices = new Point3d[19];
            double u = 0;
            double v = 0;
            double w = 0;
            u = Vsize.X / 2;
            v = Vsize.Y / 2;
            w = Vsize.Z / 2;

            Vertices[00] = VoxelPoint;
            Vertices[01] = VoxelPoint + new Vector3d(+u, 0, 0);
            Vertices[02] = VoxelPoint + new Vector3d(0, +v, 0);
            Vertices[03] = VoxelPoint + new Vector3d(0, 0, +w);
            Vertices[04] = VoxelPoint + new Vector3d(+u, +v, 0);
            Vertices[05] = VoxelPoint + new Vector3d(0, +v, +w);
            Vertices[06] = VoxelPoint + new Vector3d(+u, 0, +w);
            Vertices[07] = VoxelPoint - new Vector3d(+u, 0, 0);
            Vertices[08] = VoxelPoint - new Vector3d(0, +v, 0);
            Vertices[09] = VoxelPoint - new Vector3d(0, 0, +w);
            Vertices[10] = VoxelPoint - new Vector3d(+u, +v, 0);
            Vertices[11] = VoxelPoint - new Vector3d(0, +v, +w);
            Vertices[12] = VoxelPoint - new Vector3d(+u, 0, +w);
            Vertices[13] = VoxelPoint + new Vector3d(+u, -v, 0);
            Vertices[14] = VoxelPoint + new Vector3d(-u, +v, 0);
            Vertices[15] = VoxelPoint + new Vector3d(+u, 0, -w);
            Vertices[16] = VoxelPoint + new Vector3d(-u, 0, +w);
            Vertices[17] = VoxelPoint + new Vector3d(0, +v, -w);
            Vertices[18] = VoxelPoint + new Vector3d(0, -v, +w);

            dynamic Faces = IntersectionTarget.Faces;
            Faces.AddFace(0, 1, 4, 2);
            //xy
            Faces.AddFace(0, 2, 14, 7);
            //xy
            Faces.AddFace(0, 7, 10, 8);
            //xy
            Faces.AddFace(0, 8, 13, 1);
            //xy

            Faces.AddFace(0, 1, 6, 3);
            //xz
            Faces.AddFace(0, 3, 16, 7);
            //xZ
            Faces.AddFace(0, 7, 12, 9);
            //xz
            Faces.AddFace(0, 9, 15, 1);
            //xz

            Faces.AddFace(0, 2, 5, 3);
            //yZ
            Faces.AddFace(0, 3, 18, 8);
            //yz
            Faces.AddFace(0, 8, 11, 9);
            //yz
            Faces.AddFace(0, 9, 17, 2);
            //yz

            IntersectionTarget.Vertices.AddVertices(Vertices);
            IntersectionTarget.Faces.AddFaces(Faces);

            return IntersectionTarget;
        }
        public Mesh LineTarget6_Connected(Point3d VoxelPoint, Vector3d Vsize)
        {
            Mesh IntersectionTarget = new Mesh();
            Point3d[] Vertices = new Point3d[9];
            double u = 0;
            double v = 0;
            double w = 0;
            u = Vsize.X / 2;
            v = Vsize.Y / 2;
            w = Vsize.Z / 2;

            Vertices[0] = VoxelPoint;
            Vertices[1] = VoxelPoint - new Vector3d(+u, +v, +w);
            Vertices[2] = VoxelPoint - new Vector3d(-u, +v, +w);
            Vertices[3] = VoxelPoint - new Vector3d(-u, -v, +w);
            Vertices[4] = VoxelPoint - new Vector3d(+u, -v, +w);
            Vertices[5] = VoxelPoint + new Vector3d(-u, -v, +w);
            Vertices[6] = VoxelPoint + new Vector3d(+u, -v, +w);
            Vertices[7] = VoxelPoint + new Vector3d(+u, +v, +w);
            Vertices[8] = VoxelPoint + new Vector3d(-u, +v, +w);

            dynamic Faces = IntersectionTarget.Faces;

            Faces.AddFace(0, 1, 2);
            Faces.AddFace(0, 2, 6);
            Faces.AddFace(0, 6, 5);
            Faces.AddFace(0, 5, 1);
            Faces.AddFace(0, 4, 1);
            Faces.AddFace(0, 5, 8);
            Faces.AddFace(0, 8, 7);
            Faces.AddFace(0, 7, 3);
            Faces.AddFace(0, 3, 4);
            Faces.AddFace(0, 4, 8);
            Faces.AddFace(0, 7, 6);
            Faces.AddFace(0, 2, 3);

            IntersectionTarget.Vertices.AddVertices(Vertices);
            IntersectionTarget.Faces.AddFaces(Faces);

            return IntersectionTarget;
        }
        public List<Point3d> NearPoints(List<Point3d> x, Curve y, double Distance)
        {
            List<Point3d> AllPoints = x;
            double t = 0;
            List<Point3d> RelevantPoints = AllPoints.FindAll(lambda => y.ClosestPoint(lambda,out t, Distance));
            return RelevantPoints;
        }
        public object NearPoints(List<Point3d> x, Mesh y, double Distance)
        {
            List<Point3d> AllPoints = x;
            Point3d MP = default(Point3d);
            List<Point3d> RelevantPoints = AllPoints.FindAll(lambda => !((y.ClosestPoint(lambda, out MP, Distance))==-1));
            return RelevantPoints;
        }
    }
}
 
