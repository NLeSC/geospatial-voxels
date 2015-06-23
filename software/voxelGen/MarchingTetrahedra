//Code written by Pirouz Nourian, researcher at TU Delft 3D geoinformation  in 2014, based on the algorithm described by Paul Bourke; released under a Creative Commons license. http://creativecommons.org/licenses/by/4.0/
//http://paulbourke.net/geometry/polygonise/
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino;
using Rhino.Geometry;

namespace MarchingTetrahedrons_CSharp
{
    public class Class1
    {
        public Mesh Field_to_IsoSurface(double ts, List<Point3d> p, List<double> m, int x, int y, int z)
        {
            List<double[]> GMs = null;
            List<Box> GBoxes = Grid_to_Boxes(Plane.WorldXY, p, m, x, y, z, ref GMs);
            //Print(GBoxes.count)
            List<Tetrahedron> Tetrahedra = new List<Tetrahedron>();
            for (int k = 0; k <= GBoxes.Count - 1; k++)
            {
                Tetrahedra.AddRange(BoxToTetrahedra(GBoxes[k], GMs[k].ToArray(), ts));
            }
            Mesh IsoSurface = new Mesh();
            foreach (Tetrahedron TH in Tetrahedra)
            {
                string BSC = TH.BinaryStatus;
                if (BSC.Equals("0000") | BSC.Equals("1111"))
                {
                }
                else
                {
                    IsoSurface.Append(MeshATetrahderon(TH, ts));
                }
            }
            IsoSurface.UnifyNormals();
            IsoSurface.MakeDeformable();
            //IsoSurface.Weld(Math.PI)
            return IsoSurface;
        }
        public class Tetrahedron
        {
            private int[] indices = new int[4];
            private Point3d[] vertexPoints = new Point3d[4];
            private Point3d[] BoxVertices = new Point3d[4];
            private double[] Measures = new double[4];
            private string Status = null;
            public Tetrahedron(Box ABox, int[] Vertex_Indices, double[] BoxMeasures, double isoValue)
            {
                if (BoxMeasures.Length != 8)
                    return;
                if (Vertex_Indices.Length != 4)
                    return;
                BoxVertices = ABox.GetCorners();
                if (BoxVertices.Length != 8)
                    return;
                indices = Vertex_Indices;
                for (int i = 0; i <= 3; i++)
                {
                    vertexPoints[i] = BoxVertices[Vertex_Indices[i]];
                    Measures[i] = BoxMeasures[Vertex_Indices[i]];
                    if (Measures[i] > isoValue)
                    {
                        Status = "1" + Status;
                    }
                    else
                    {
                        Status = "0" + Status;
                    }
                }
            }
            public double[] VertexMeasures
            {
                get { return Measures; }
                set { Measures = value; }
            }
            public Point3d[] Vertices
            {
                get { return this.vertexPoints; }
            }
            public string Description
            {
                get { return string.Join(",", indices); }
            }
            public string BinaryStatus
            {
                get { return Status; }
            }
            public Mesh toMesh()
            {
                Mesh TetraMesh = new Mesh();
                TetraMesh.Vertices.AddVertices(this.vertexPoints);
                TetraMesh.Faces.AddFace(new MeshFace(0, 2, 1));TetraMesh.Faces.AddFace(new MeshFace(1, 2, 3));
                TetraMesh.Faces.AddFace(new MeshFace(0, 1, 3));TetraMesh.Faces.AddFace(new MeshFace(0, 3, 2));
                return TetraMesh;
            }
        }
        public List<Tetrahedron> BoxToTetrahedra(Box ABox, double[] Measures, double iso)
        {
            if (Measures.Length != 8)
                return null;
            List<Tetrahedron> tetrahedra = new List<Tetrahedron>();
            tetrahedra.Add(new Tetrahedron(ABox, new int[] { 0, 2, 3, 7 }, Measures, iso));tetrahedra.Add(new Tetrahedron(ABox, new int[] { 0, 2, 6, 7 }, Measures, iso));
            tetrahedra.Add(new Tetrahedron(ABox, new int[] { 0, 4, 6, 7 }, Measures, iso));tetrahedra.Add(new Tetrahedron(ABox, new int[] { 0, 6, 1, 2 }, Measures, iso));
            tetrahedra.Add(new Tetrahedron(ABox, new int[] { 0, 6, 1, 4 }, Measures, iso));tetrahedra.Add(new Tetrahedron(ABox, new int[] { 5, 6, 1, 4 }, Measures, iso));
            return tetrahedra;
        }
        public List<Box> Grid_to_Boxes(Plane BasePlane, List<Point3d> GridPoints, List<double> GridMeasures, int XCount, int YCount, int ZCount, ref List<double[]> BoxMs)
        {
            if (GridPoints.Count != GridMeasures.Count)
                return null;
            List<Box> GBoxes = new List<Box>();
            List<double[]> GMs = new List<double[]>();
            int i = 0;
            int j = 0;
            int k = 0;
            for (k = 0; k <= ZCount - 2; k++)
            {
                for (j = 0; j <= YCount - 2; j++)
                {
                    for (i = 0; i <= XCount - 2; i++)
                    {
                        string code = (k + j + i).ToString();
                        //print(code)
                        List<Point3d> points = new List<Point3d>();
                        double[] GMeasures = new double[8];
                        int n0 = 0;int n1 = 0;int n2 = 0;int n3 = 0;int n4 = 0;int n5 = 0;int n6 = 0;int n7 = 0;
                        n0 = i + XCount * j + XCount * YCount * k;
                        n1 = n0 + 1;
                        n2 = n0 + XCount;
                        n3 = n1 + XCount;
                        n4 = n0 + XCount * YCount;
                        n5 = n4 + 1;
                        n6 = n4 + XCount;
                        n7 = n6 + 1;
                        points.Add(GridPoints[n0]);points.Add(GridPoints[n1]);points.Add(GridPoints[n2]);points.Add(GridPoints[n3]);
                        points.Add(GridPoints[n4]);points.Add(GridPoints[n5]);points.Add(GridPoints[n6]);points.Add(GridPoints[n7]);
                        Box GB = new Box(BasePlane, points);
                        GMeasures = new double[]
              {
              GridMeasures[n0],GridMeasures[n1],GridMeasures[n3],GridMeasures[n2],GridMeasures[n4],GridMeasures[n5],GridMeasures[n7],GridMeasures[n6]
              };
                        //the order does not seem reasonable but it is. This is because of the way the box is formed out of vertices
                        GMs.Add(GMeasures);
                        GBoxes.Add(GB);
                    }
                }
            }
            BoxMs = GMs;
            return GBoxes;
        }
        public Mesh MeshATetrahderon(Tetrahedron TH, double ts)
        {
            if (TH.VertexMeasures.Length != 4)
                return null;
            object BSC = TH.BinaryStatus;
            Mesh resMesh = new Mesh();
            Mesh IsoMesh = new Mesh();
            List<Point3d> Vertices = new List<Point3d>();
            MeshFace Face = default(MeshFace);
            if (BSC.Equals("0000") | BSC.Equals("1111"))
            {
                return null;
            }
            else if (BSC.Equals("0001") | BSC.Equals("1110"))
            {
                Vertices.Add(InterPoint(0, 1, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(0, 2, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(0, 3, ts, TH.Vertices, TH.VertexMeasures));
                Face = new MeshFace(0, 1, 2);
            }
            else if (BSC.Equals("0010") | BSC.Equals("1101"))
            {
                Vertices.Add(InterPoint(1, 0, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(1, 2, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(1, 3, ts, TH.Vertices, TH.VertexMeasures));
                Face = new MeshFace(2, 1, 0);
            }
            else if (BSC.Equals("0100") | BSC.Equals("1011"))
            {
                Vertices.Add(InterPoint(2, 0, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(2, 1, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(2, 3, ts, TH.Vertices, TH.VertexMeasures));
                Face = new MeshFace(2, 1, 0);
            }
            else if (BSC.Equals("1000") | BSC.Equals("0111"))
            {
                Vertices.Add(InterPoint(3, 0, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(3, 1, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(3, 2, ts, TH.Vertices, TH.VertexMeasures));
                Face = new MeshFace(2, 1, 0);
            }
            else if (BSC.Equals("0011") | BSC.Equals("1100"))
            {
                Vertices.Add(InterPoint(0, 2, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(0, 3, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(1, 2, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(1, 3, ts, TH.Vertices, TH.VertexMeasures));
                Face = new MeshFace(0, 2, 3, 1);
            }
            else if (BSC.Equals("0101") | BSC.Equals("1010"))
            {
                Vertices.Add(InterPoint(0, 1, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(0, 3, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(2, 1, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(2, 3, ts, TH.Vertices, TH.VertexMeasures));
                Face = new MeshFace(0, 2, 3, 1);
            }
            else if (BSC.Equals("0110") | BSC.Equals("1001"))
            {
                Vertices.Add(InterPoint(1, 0, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(1, 3, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(2, 0, ts, TH.Vertices, TH.VertexMeasures));
                Vertices.Add(InterPoint(2, 3, ts, TH.Vertices, TH.VertexMeasures));
                Face = new MeshFace(0, 2, 3, 1);
            }
            resMesh.Vertices.AddVertices(Vertices);
            resMesh.Faces.AddFace(Face);
            return resMesh;
        }
        public Point3d InterPoint(int p1, int p2, double iso, Point3d[] P, double[] M)
        {
            Point3d IntP = default(Point3d);IntP = P[p1] + ((P[p2] - P[p1]) * ((iso - M[p1]) / (M[p2] - M[p1])));
            return IntP;
        }
    }
}
