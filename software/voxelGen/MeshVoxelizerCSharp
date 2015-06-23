//Algorithm developed and written by Pirouz Nourian, researcher @ TU Delft, 3D Geoinformation in 2014. 
//The work reported here has been funded by NLeSC Big Data Analytics in the Geo-Spatial
using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using Rhino;
using Rhino.Geometry;

namespace TUD_TopoVoxelizer
{
    public class VoxelizationFunctions
    {  
        public List<Tuple<int, double, double, double, byte, byte, byte>> VoxelizeMesh(Mesh M, Color MC, int M_ID, Vector3d Vsize, int CO, ref PointCloud VPC)
        {
            //{ID,X,Y,Z,R,G,B}
            BoundingBox RBBox = M.GetBoundingBox(false);
            //bounding box in R3 (real numbers)
            BoundingBox ZBBox = RBBox_to_ZBBox(RBBox, Vsize);
            //boudning box in Z3 (integer numbers) (note that they are again embedded in R3 for visualization)

            List<Point3d> gridpoints = BBoxToVoxels(ZBBox, Vsize);

            List<Point3d> Mesh_NearPoints = NearPoints(gridpoints, M, Vsize.Length / 2);
            //A = Mesh_NearPoints
            List<Tuple<int, double, double, double, byte, byte, byte>> Voxels = new List<Tuple<int, double, double, double, byte, byte, byte>>();
            //{ID,X,Y,Z,R,G,B}

            Line[] Intarget = null;
            //int[] Facets = null;
            // for each potential voxel check if it has to be included in the raster
            foreach (Point3d gp in Mesh_NearPoints)
            {
                if (CO == 26)
                {
                    Intarget = MeshTarget26_Connected(gp, Vsize);
                }
                else if (CO == 6)
                {
                    Intarget = MeshTarget06_Connected(gp, Vsize);
                }
                else
                {
                    RhinoApp.WriteLine("connectivity target undefined!");
                }
                //Line[] Intersections = Array.FindAll(Intarget, lambda => MeshLineIntersect(M, lambda));
                if (Array.Exists(Intarget, lambda => MeshLineIntersect(M, lambda)))
                {
                    Voxels.Add(Tuple.Create(M_ID, gp.X, gp.Y, gp.Z, MC.R, MC.G, MC.B));
                    VPC.Add(gp, MC);
                }
            }
            return Voxels;
        }
        public Mesh[] MeshToTriangles(Mesh M) //not necessary perhaps; just in case... a utility function
        {
            M.Faces.ConvertQuadsToTriangles();
            Mesh[] Meshes = null;
            M.Unweld(0, false);
            Meshes = M.ExplodeAtUnweldedEdges();
            return Meshes;
        }
        public BoundingBox RBBox_to_ZBBox(BoundingBox RBbox, Vector3d VSize) // converts a bounding box given in R3 to a bounding box in Z3 that is bigger than or the same size as the one given in R3
        {
            Point3d ZMinP = MinBoundRP_to_ZP(RBbox.Min, VSize);
            Point3d ZMaxP = MaxBoundRP_to_ZP(RBbox.Max, VSize);
            BoundingBox ZBBox = new BoundingBox(ZMinP, ZMaxP);
            return ZBBox;
        }
        public Point3d RPtoZP(Point3d RPoint, Vector3d VSize) // finds the closest Z3 point to an R3 point
        {
            Point3d ZPOint = new Point3d((Math.Floor(RPoint.X / VSize.X) + 0.5) * VSize.X, (Math.Floor(RPoint.Y / VSize.Y) + 0.5) * VSize.Y, (Math.Floor(RPoint.Z) / VSize.Z + 0.5) * VSize.Z);
            return ZPOint;
        }
        private Point3d MaxBoundRP_to_ZP(Point3d RPoint, Vector3d VSize) //snaps a boundingbox Maximum Corner point embedded in R3 to a point in Z3 ensuring a bounding box in Z3 that is bigger than or the same size as the bounding box in R3
        {
            double x = 0; double y = 0; double z = 0;
            double u = 0; double v = 0; double w = 0;
            double ZP_x = 0;
            double ZP_y = 0;
            double ZP_z = 0;
            x = RPoint.X;
            y = RPoint.Y;
            z = RPoint.Z;
            u = VSize.X;
            v = VSize.Y;
            w = VSize.Z;
            ZP_x = Math.Ceiling(x / u) * u;
            ZP_y = Math.Ceiling(y / v) * v;
            ZP_z = Math.Ceiling(z / w) * w;
            Point3d ZPoint = new Point3d(ZP_x, ZP_y, ZP_z);
            return ZPoint;
        }
        private Point3d MinBoundRP_to_ZP(Point3d RPoint, Vector3d VSize) //snaps a boundingbox Minimum Corner point embedded in R3 to a point in Z3 ensuring a bounding box in Z3 that is bigger than or the same size as the bounding box in R3
        {
            double x = 0; double y = 0; double z = 0;
            double u = 0; double v = 0; double w = 0;
            double ZP_x = 0;
            double ZP_y = 0;
            double ZP_z = 0;
            x = RPoint.X;
            y = RPoint.Y;
            z = RPoint.Z;
            u = VSize.X;
            v = VSize.Y;
            w = VSize.Z;
            ZP_x = Math.Floor(x / u) * u;
            ZP_y = Math.Floor(y / v) * v;
            ZP_z = Math.Floor(z / w) * w;
            Point3d ZPoint = new Point3d(ZP_x, ZP_y, ZP_z);
            return ZPoint;
        }
        public List<Point3d> BBoxToVoxels(BoundingBox ZBBox, Vector3d Vsize) //Creates a List of voxels that correspond to a boundingbox described in Z3 space
        {
            int i = 0;
            int j = 0;
            int k = 0;
            int Imax = 0;
            int Jmax = 0;
            int Kmax = 0;
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
        public Line[] MeshTarget26_Connected(Point3d VoxelPoint, Vector3d Vsize) ////makes a target for intersection with a mesh that ensures 26-connected results as proven by Samuli Laine 2013, NVIDIA Research: this is the set of a cube's diagons
        {
            //3D crosshair target for 26-connected results
            Line[] IntersectionTarget = new Line[3];
            Point3d[] Vertices = new Point3d[6];
            double u = 0;
            double v = 0;
            double w = 0;
            u = Vsize.X / 2;
            v = Vsize.Y / 2;
            w = Vsize.Z / 2;
            Vertices[0] = VoxelPoint + new Vector3d(+u, 0, 0);
            Vertices[1] = VoxelPoint + new Vector3d(0, +v, 0);
            Vertices[2] = VoxelPoint + new Vector3d(0, 0, +w);
            Vertices[3] = VoxelPoint + new Vector3d(-u, 0, 0);
            Vertices[4] = VoxelPoint + new Vector3d(0, -v, 0);
            Vertices[5] = VoxelPoint + new Vector3d(0, 0, -w);
            for (int i = 0; i <= 2; i++)
            {
                IntersectionTarget[i] = new Line(Vertices[i], Vertices[i + 3]);
            }
            return IntersectionTarget;
        }
        public Line[] MeshTarget06_Connected(Point3d VP, Vector3d Vsize) //makes a target for intersection with a mesh that ensures 06-connected results as proven by Samuli Laine 2013, NVIDIA Research: this is the outline wireframe of a cube
        {
            //cube outline target for 6-connected results
            Line[] IT = new Line[12];
            Point3d[] VX = new Point3d[8];
            double u = 0;
            double v = 0;
            double w = 0;
            u = Vsize.X / 2;
            v = Vsize.Y / 2;
            w = Vsize.Z / 2;

            VX[0] = VP + new Vector3d(+u, +v, +w);
            VX[1] = VP + new Vector3d(-u, +v, +w);
            VX[2] = VP + new Vector3d(-u, -v, +w);
            VX[3] = VP + new Vector3d(+u, -v, +w);
            VX[4] = VP + new Vector3d(-u, -v, -w);
            VX[5] = VP + new Vector3d(+u, -v, -w);
            VX[6] = VP + new Vector3d(+u, +v, -w);
            VX[7] = VP + new Vector3d(-u, +v, -w);
            //cube edges
            IT[00] = new Line(VX[0], VX[1]);
            IT[01] = new Line(VX[1], VX[2]);
            IT[02] = new Line(VX[2], VX[3]);
            IT[03] = new Line(VX[3], VX[0]);
            IT[04] = new Line(VX[0], VX[6]);
            IT[05] = new Line(VX[6], VX[5]);
            IT[06] = new Line(VX[5], VX[4]);
            IT[07] = new Line(VX[4], VX[7]);
            IT[08] = new Line(VX[5], VX[3]);
            IT[09] = new Line(VX[4], VX[2]);
            IT[10] = new Line(VX[1], VX[7]);
            IT[11] = new Line(VX[6], VX[7]);
            return IT;
        }
        public List<Point3d> NearPoints(List<Point3d> x, Mesh y, double Distance) //finds the voxel centrepoints that could be relevant as to their distance towards the mesh, which is to be voxelized.
        {
            List<Point3d> AllPoints = x;
            Point3d MP = default(Point3d);
            List<Point3d> RelevantPoints = AllPoints.FindAll(lambda => y.ClosestPoint(lambda, out MP, Distance) != -1);
            return RelevantPoints;
        }
        public List<Point3d> NearPoints(List<Point3d> x, Curve y, double Distance) //finds the voxel centrepoints that could be relevant as to their distance towards the curve, which is to be voxelized.
        {
            List<Point3d> AllPoints = x;
            double t = 0;
            List<Point3d> RelevantPoints = AllPoints.FindAll(lambda => y.ClosestPoint(lambda, out t, Distance));
            return RelevantPoints;
        }
        public bool MeshLineIntersect(Rhino.Geometry.Mesh x, Rhino.Geometry.Line L) //this is now for reading triangles from Rhino Mesh objects; can be eventually replaced with something, which reads OBJ mesh objects
        {
            //List<Point3d[]> TVs = new List<Point3d[]>();
            //List of Triangle Vertices
            Point3d[] TV = new Point3d[3];
            //Triangle Vertices
            Point3f av = default(Point3f);
            Point3f bv = default(Point3f);
            Point3f cv = default(Point3f);
            Point3f dv = default(Point3f);

            //Each face As MeshFace In x.Faces
            bool does = false;
            for (int k = 0; k <= x.Faces.Count - 1; k++)
            {
                x.Faces.GetFaceVertices(k, out av, out bv, out cv, out dv);
                TV = new Point3d[] { new Point3d(av), new Point3d(bv), new Point3d(cv) };
                //TVs.Add(TV);
                if (TriangleLineIntersect(TV, x.FaceNormals[k], L))
                {
                    does = true;
                }
            }
            return does;//TVs.Exists(Lambda => TriangleLineIntersect(Lambda, L));
        }
        public bool TriangleLineIntersect(Rhino.Geometry.Point3d[] Vx, Rhino.Geometry.Line L)
        {
            if (Vx.Length != 3)
            {
                throw new Exception("Triangle is not valid!");
            }
            Point3d O = Vx[0];
            Vector3d U = Vx[1] - O;
            Vector3d V = Vx[2] - O;
            Vector3d N = Vector3d.CrossProduct(U, V);
            //plane normal

            Point3d PS = L.From;//start point of the line
            Point3d PE = L.To;  //end point of the line

            double Nomin = ((O - PS) * N);
            //operator * for dot product
            double Denom = N * (PE - PS);

            // only if the line is not paralell to the plane containing triangle T
            if (Denom != 0)
            {
                double alpha = Nomin / Denom;
                // parameter along the line where it intersects the plane in question, only if not paralell to the plane
                Point3d P = PS + alpha * (PE - PS);
                //L.PointAt(alpha) '

                Vector3d W = P - O;

                double UU = U * U;
                double VV = V * V;
                double UV = U * V;
                double WU = W * U;
                double WV = W * V;

                double STDenom = Math.Pow(UV, 2) - UU * VV;

                double s = (UV * WV - VV * WU) / STDenom;
                double t = (UV * WU - UU * WV) / STDenom;

                Point3d Point = O + s * U + t * V;

                if (s >= 0 & t >= 0 & s + t <= 1)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        public bool TriangleLineIntersect(Rhino.Geometry.Point3d[] Vx, Rhino.Geometry.Vector3d Normal, Rhino.Geometry.Line L)
        {
            if (Vx.Length != 3)
            {
                throw new Exception("Triangle is not valid!");
            }
            Point3d O = Vx[0];
            Vector3d U = Vx[1] - O;
            Vector3d V = Vx[2] - O;
            Vector3d N = Normal;//Vector3d.CrossProduct(U, V);
            //plane normal

            Point3d PS = L.From;//start point of the line
            Point3d PE = L.To;  //end point of the line

            double Nomin = ((O - PS) * N);
            //operator * for dot product
            double Denom = N * (PE - PS);

            // only if the line is not paralell to the plane containing triangle T
            if (Denom != 0)
            {
                double alpha = Nomin / Denom;
                // parameter along the line where it intersects the plane in question, only if not paralell to the plane
                Point3d P = PS + alpha * (PE - PS);
                //L.PointAt(alpha) '

                Vector3d W = P - O;

                double UU = U * U;
                double VV = V * V;
                double UV = U * V;
                double WU = W * U;
                double WV = W * V;

                double STDenom = Math.Pow(UV, 2) - UU * VV;

                double s = (UV * WV - VV * WU) / STDenom;
                double t = (UV * WU - UU * WV) / STDenom;

                Point3d Point = O + s * U + t * V;

                if (s >= 0 & t >= 0 & s + t <= 1)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        public void VoxelTextFileWriter(int CO, List<Mesh> ObjM, List<Color> ObjC, string FPath, double Sx, double Sy, double Sz, ref object A)
        {
            List<Mesh> Meshes = ObjM;
            List<Color> Colors = ObjC;

            Vector3d VSize = new Vector3d(Sx, Sy, Sz);

            double Distance = VSize.Length / 2;

            List<Tuple<int, double, double, double, byte, byte, byte>> VT = new List<Tuple<int, double, double, double, byte, byte, byte>>();
            //Voxel Tuples'List(Of String())'
            PointCloud VPC = new PointCloud();
            //VoxelPointCloud
            for (int i = 0; i <= ObjM.Count - 1; i++)
            {
                VT.AddRange(VoxelizeMesh(Meshes[i], Colors[i], i, VSize, CO, ref VPC));
            }

            System.IO.StreamWriter SW = File.CreateText(FPath);
            string result = "";
            //k As Integer=0 To VT.count - 1'String() In VT'
            foreach (Tuple<int, double, double, double, byte, byte, byte> vtp in VT)
            {
                result = vtp.Item1 + "," + vtp.Item2 + "," + vtp.Item3 + "," + vtp.Item4 + "," + vtp.Item5 + "," + vtp.Item6 + "," + vtp.Item7;
                //String.Join(",", vtp) 'vtp.toString'
                SW.WriteLine(result);
            }
            SW.Close();

            //PCloud = null;
            //PCloud = VPC;

        }
    }
    public class independentGeoMethods
    {
        public double Dist_Point_Segment(Point3d P, Line S)
        {
            Point3d P0 = S.From;
            //start point of segemnt
            Point3d P1 = S.To;
            //end poiunt of segment
            Vector3d V = S.To - P0;
            //EndPoint-StartPoint
            Vector3d W = P - P0;
            //the vector from point to start of the segment
            double c1 = W * V;
            // * is used for dot product
            if (c1 <= 0)
                return P.DistanceTo(P0);
            //Eucleadian distance of the point to the start of the line
            double c2 = V * V;
            if (c2 <= c1)
                return P.DistanceTo(P1);
            //Eucleadian distance of the point to the end of the line
            double b = c1 / c2;
            Point3d Pb = P0 + b * V;
            //point on the line segment closest to P
            return P.DistanceTo(Pb);
        }

        public double Dist_point_Triangle(Point3d[] Vx, Point3d P)
        {
            Point3d O = Vx[0];
            Vector3d U = Vx[1] - O;
            Vector3d V = Vx[2] - O;

            Vector3d W = P - O;

            double UU = U * U;
            double VV = V * V;
            double UV = U * V;
            double WU = W * U;
            double WV = W * V;

            double det = Math.Abs(Math.Pow(UV, 2) - UU * VV);
            // determinant will be >0 only if the two edges are not colinear (i.e. the triangle is NOT degenerate)
            if (det > 0)
            {
                double s = -(UV * WV - VV * WU);
                /// det
                double t = -(UV * WU - UU * WV);
                /// det

                //left half of the line(V1,V2)
                if (s + t <= det)
                {
                    if (s < 0)
                    {
                        if (t < 0)
                        {
                            //Region4, distance to V0
                            return Vx[0].DistanceTo(P);            //P.DistanceTo(Vx(0))
                        }
                        else
                        {
                            //Region3, distance to line(V0,V2)
                            return Dist_Point_Segment(P, new Line(Vx[0], Vx[2]));
                        }
                    }
                    else if (t < 0)
                    {
                        //Region5, distance to line(V0,V1)

                        return Dist_Point_Segment(P, new Line(Vx[0], Vx[1]));
                    }
                    else
                    {
                        //Region0, distance to point(s,t)

                        s = s / det;
                        t = t / det;
                        Point3d Vertex = O + s * U + t * V;
                        return Vertex.DistanceTo(P);

                    }
                    // right half of the line(V1,V2)
                }
                else
                {
                    if (s < 0)
                    {
                        //Region2, distance to V2

                        return Vx[2].DistanceTo(P);
                        //P.DistanceTo(Vx(2))
                    }
                    else if (t < 0)
                    {
                        //Region6, distance to V1

                        return Vx[1].DistanceTo(P);
                        //P.DistanceTo(Vx(1))
                    }
                    else
                    {
                        //Region1, distance to line(V1,V2)

                        return Dist_Point_Segment(P, new Line(Vx[1], Vx[2]));
                    }
                }
            }
            else
            {
                return -1;//det=0, the triangle is degenerate!!!
                throw new Exception("there is a degenerate triangle!");
            }
        }

            }

    }
 
