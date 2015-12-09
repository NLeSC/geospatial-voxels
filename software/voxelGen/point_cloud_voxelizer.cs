//Algorithm developed and written by Pirouz Nourian, researcher @ TU Delft, 3D Geoinformation in 2014. 
//The work reported here has been funded by NLeSC Big Data Analytics in the Geo-Spatial
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;


namespace PointCloudVoxelizer
{
    public class Class1
    {
        private void RunScript(Guid id, double Sx, double Sy, double Sz, ref object A, ref object B)
        {
            object obj = RhinoDocument.Objects.Find(id);
            Rhino.DocObjects.PointCloudObject cloudObject = obj as Rhino.DocObjects.PointCloudObject;
            if (cloudObject == null)
            {
                return;
            }
            PointCloud cloud = cloudObject.PointCloudGeometry;
            //    var PCO = new { PC = cloud };
            //    A = PCO;
            //cloudObject.PointCloudGeometry
            //Get Bounding Box of Mesh and find the X,Y,Z extents discritizations:
            BoundingBox BB = cloud.GetBoundingBox(false);
            //M.GetBoundingBox(True)
            Vector3d VSize = new Vector3d(Sx, Sy, Sz);
            BoundingBox Zbbox = RBBox_to_ZBBox(BB, VSize);
            int CX = Convert.ToInt32(Math.Floor((Zbbox.Max.X - Zbbox.Min.X) / Sx));
            int CY = Convert.ToInt32(Math.Floor((Zbbox.Max.Y - Zbbox.Min.Y) / Sy));
            int CZ = Convert.ToInt32(Math.Floor((Zbbox.Max.Z - Zbbox.Min.Z) / Sz));
            bool[, ,] Voxels3I = new bool[CX + 1, CY + 1, CZ + 1];
            List<Point3d> VoxelsOut = new List<Point3d>();
            //FileSystem.Print(string.Format("{0},{1},{2}", CX, CY, CZ));
            foreach (Point3d point in cloud.GetPoints())
            {
                int i = Convert.ToInt32(Math.Floor((point.X - Zbbox.Min.X) / Sx));
                int j = Convert.ToInt32(Math.Floor((point.Y - Zbbox.Min.Y) / Sy));
                int k = Convert.ToInt32(Math.Floor((point.Z - Zbbox.Min.Z) / Sz));
                if (!Voxels3I[i, j, k])
                {
                    Voxels3I[i, j, k] = true;
                    Point3d VoxelOut = new Point3d(i * Sx, j * Sy, k * Sz) + VSize / 2 + Zbbox.Min;
                    VoxelsOut.Add(VoxelOut);
                }
            }

            A = VoxelsOut;
            B = Voxels_to_Boxels(VoxelsOut.ToArray(), new Vector3d(Sx, Sy, Sz));
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
            ZP_x = Math.Ceiling(x / u) * u;
            ZP_y = Math.Ceiling(y / v) * v;
            ZP_z = Math.Ceiling(z / w) * w;
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
            ZP_x = Math.Floor(x / u) * u;
            ZP_y = Math.Floor(y / v) * v;
            ZP_z = Math.Floor(z / w) * w;
            Point3d ZPoint = new Point3d(ZP_x, ZP_y, ZP_z);
            return ZPoint;
        }
        public object Voxels_to_Boxels(Point3d[] VoxelPoints, Vector3d VSize)
        {
            List<Mesh> Cubes = new List<Mesh>();
            double u = 0;            double v = 0;            double w = 0;
            u = VSize.X / 2;            v = VSize.Y / 2;            w = VSize.Z / 2;

            foreach (Point3d voxelPoint in VoxelPoints)
            {
                Point3d xp = voxelPoint;
                Point3d[] VX = new Point3d[8];
                VX[0] = xp - new Vector3d(+u, +v, +w);
                VX[1] = xp - new Vector3d(-u, +v, +w);
                VX[2] = xp - new Vector3d(-u, -v, +w);
                VX[3] = xp - new Vector3d(+u, -v, +w);
                VX[4] = xp + new Vector3d(-u, -v, +w);
                VX[5] = xp + new Vector3d(+u, -v, +w);
                VX[6] = xp + new Vector3d(+u, +v, +w);
                VX[7] = xp + new Vector3d(-u, +v, +w);
                Mesh Cube = Mesh.CreateFromBox(VX, 1, 1, 1);
                //Cube.VertexColors.CreateMonotoneMesh(Voxel.SemanticToColor(SemanticTags, SemanticColors, vxl.SemanticIndex))
                //Cube.VertexColors.CreateMonotoneMesh(Color.black)
                Cubes.Add(Cube);
            }
            return Cubes;
        }
    }
}

