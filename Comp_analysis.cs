using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using ScottPlot.Statistics;
using static WindowsFormsApp3.Form1;

namespace WindowsFormsApp3
{


    public class Comp_analysis
    {
        static double angle_RAD;
        public static double[,] histogram;
        public static void to_independent_2D(Form1 my_form)
        {
            angle_RAD = Math.Atan((2.0 * my_form.r * my_form.middle_squared_multi[0] * my_form.middle_squared_multi[1]) / (Math.Pow(my_form.middle_squared_multi[0], 2) - Math.Pow(my_form.middle_squared_multi[1], 2)));
            angle_RAD /= 2.0;
            MessageBox.Show($"{angle_RAD}");
            List<double[]> data_2D_orig = new List<double[]>();
            for (int i = 0; i < my_form.all_2n.Count; i++)
            {
                double[] temp_arr = new double[my_form.all_2n[i].Length];
                temp_arr[0] = my_form.all_2n[i][my_form.vector1];
                temp_arr[1] = my_form.all_2n[i][my_form.vector2];
                data_2D_orig.Add(temp_arr);
                double x = my_form.all_2n[i][my_form.vector1];
                double y = my_form.all_2n[i][my_form.vector2];
                my_form.all_2n[i][my_form.vector1] = x * Math.Cos(angle_RAD) + y * Math.Sin(angle_RAD);
                my_form.all_2n[i][my_form.vector2] = -1.0 * x * Math.Sin(angle_RAD) + y * Math.Cos(angle_RAD);
                //my_form.all_2n[i][my_form.vector1] = x * Math.Cos(angle_RAD) - y * Math.Sin(angle_RAD);
                //my_form.all_2n[i][my_form.vector2] = x * Math.Sin(angle_RAD) + y * Math.Cos(angle_RAD);
                if (i == 0)
                {
                    my_form.max_n2[0] = my_form.all_2n[i][my_form.vector1];
                    my_form.max_n2[1] = my_form.all_2n[i][my_form.vector2];
                    my_form.min_n2[0] = my_form.all_2n[i][my_form.vector1];
                    my_form.min_n2[1] = my_form.all_2n[i][my_form.vector2];
                }
                else
                {
                    if (my_form.max_n2[0] < my_form.all_2n[i][my_form.vector1])
                    {
                        my_form.max_n2[0] = my_form.all_2n[i][my_form.vector1];
                    }
                    if (my_form.max_n2[1] < my_form.all_2n[i][my_form.vector2])
                    {
                        my_form.max_n2[1] = my_form.all_2n[i][my_form.vector2];
                    }
                    if (my_form.min_n2[0] > my_form.all_2n[i][my_form.vector1])
                    {
                        my_form.min_n2[0] = my_form.all_2n[i][my_form.vector1];
                    }
                    if (my_form.min_n2[1] > my_form.all_2n[i][my_form.vector2])
                    {
                        my_form.min_n2[1] = my_form.all_2n[i][my_form.vector2];
                    }
                }

            }
            ;
            my_form.processing();
            histogram = new double[my_form.classes, my_form.classes];
            List<double> for_max = new List<double>();
            for (int i = 0; i < my_form.classes; i++)
            {
                for (int k = 0; k < my_form.classes; k++)
                {
                    histogram[i, k] = my_form.freq_for_histogram[i, k];
                    for_max.Add(histogram[i, k]);
                }
            }
            double max_frq = for_max.Max();
            double min_frq = for_max.Min();

            for (int i = 0; i < my_form.classes; i++)
            {
                for (int k = 0; k < my_form.classes; k++)
                {
                    histogram[i, k] = histogram[i, k] / max_frq;
                    
                }
            }


            // N == 100
            histo_X = new double[N];
            histo_Y = new double[N];
            stepX = (my_form.max_n2[0] - my_form.min_n2[0]) / (N - 1.0);
            stepY = (my_form.max_n2[1] - my_form.min_n2[1]) / (N - 1.0);
            int index = 0;
            for (double i = my_form.min_n2[0]; i <= my_form.max_n2[0]; i += stepX)
            {
                histo_X[index] = i;
                index++;
            }
            index = 0;
            for (double i = my_form.min_n2[1]; i <= my_form.max_n2[1]; i += stepY)
            {
                histo_Y[index] = i;
                index++;
            }
            pointS_Hist = new point_hist[N, N];//points
            double h_x = (my_form.max_n2[0] - my_form.min_n2[0]) / (my_form.classes);
            double h_y = (my_form.max_n2[1] - my_form.min_n2[1]) / (my_form.classes);
           
            for (int a = 0; a < N; a++)
            {
                for (int b = 0; b < N; b++)
                {


                    ////////////////////////////////
                    int ind_i = 0;
                    int ind_j = 0;
                   // for (int i = 0; i < my_form.all_2n.Count; i++)
                   // {
                        int k = 1;
                        for (; k <= my_form.classes; k++)
                        {
                            if (Math.Round(histo_X[a] - ((k - 1) * h_x + my_form.min_n2[0]), 3) >= 0 && Math.Round(histo_X[a] - (k * h_x + my_form.min_n2[0]), 3) <= 0)
                            {
                                ind_i = k - 1;
                                break;
                            }
                        }
                        k = 1;

                        for (; k <= my_form.classes; k++)
                        {
                            if (Math.Round(histo_Y[b] - ((k - 1) * h_y + my_form.min_n2[1]), 3) >= 0 && Math.Round(histo_Y[b] - (k * h_y + my_form.min_n2[1]), 3) <= 0)
                            {
                                ind_j = k - 1;
                                break;
                            }
                        }
                    //}
                    //////////////////////////////////////


                    double frq = histogram[ind_i, ind_j];
                    if(frq > 0.000000000000000)
                    {
                        ;
                    }
                    int red = Convert.ToByte((255.0 * (1.0 - frq)));
                    int green = Convert.ToByte(255.0 - frq * 15);
                    int blue = Convert.ToByte((255.0 * (1.0 - frq)));
                    Color ccccc = Color.FromArgb(255, red, green, blue);
                    //Color ccccc = Color.FromArgb(5, 250, 1);
                    //Color ccccc = Color.FromArgb((int)(10.0 * histogram[ind_i, ind_j] * 255), (int)(10.0 * histogram[ind_i, ind_j] * 255), (int)(10.0 * histogram[ind_i, ind_j] * 255));
                  
                    pointS_Hist[a,b] = new point_hist(ccccc, histo_X[a], histo_Y[b]);
                }
            }
            ;
        }
        public class point_hist
        {
            public Color colorr;
            public double x;
            public double y;
            public point_hist(Color colorr, double x, double y)
            {
                this.colorr = colorr;
                this.x = x;
                this.y = y;
            }
        }
        static int N = 100;
        static double stepX;
        static double stepY;
        static double[] histo_X;
        static double[] histo_Y;
        static point_hist[,] pointS_Hist;
        public static void to_orig_2D(Form1 my_form)
        {
            for (int i = 0; i < my_form.all_2n.Count; i++)
            {
                double x = my_form.all_2n[i][my_form.vector1];
                double y = my_form.all_2n[i][my_form.vector2];
                //my_form.all_2n[i][my_form.vector1] = x * Math.Cos(angle_RAD) + y * Math.Sin(angle_RAD);
                //my_form.all_2n[i][my_form.vector2] = -1.0 * x * Math.Sin(angle_RAD) + y * Math.Cos(angle_RAD);
                my_form.all_2n[i][my_form.vector1] = x * Math.Cos(angle_RAD) - y * Math.Sin(angle_RAD);
                my_form.all_2n[i][my_form.vector2] = x * Math.Sin(angle_RAD) + y * Math.Cos(angle_RAD);
                if (i == 0)
                {
                    my_form.max_n2[0] = my_form.all_2n[i][my_form.vector1];
                    my_form.max_n2[1] = my_form.all_2n[i][my_form.vector2];
                    my_form.min_n2[0] = my_form.all_2n[i][my_form.vector1];
                    my_form.min_n2[1] = my_form.all_2n[i][my_form.vector2];
                }
                else
                {
                    if (my_form.max_n2[0] < my_form.all_2n[i][my_form.vector1])
                    {
                        my_form.max_n2[0] = my_form.all_2n[i][my_form.vector1];
                    }
                    if (my_form.max_n2[1] < my_form.all_2n[i][my_form.vector2])
                    {
                        my_form.max_n2[1] = my_form.all_2n[i][my_form.vector2];
                    }
                    if (my_form.min_n2[0] > my_form.all_2n[i][my_form.vector1])
                    {
                        my_form.min_n2[0] = my_form.all_2n[i][my_form.vector1];
                    }
                    if (my_form.min_n2[1] > my_form.all_2n[i][my_form.vector2])
                    {
                        my_form.min_n2[1] = my_form.all_2n[i][my_form.vector2];
                    }
                }

            }
            //for (int i = 0; i < N; i++)
            //{
            //    double x = histo_X[i];
            //    double y = histo_Y[i];

            //    histo_X[i] = x * Math.Cos(angle_RAD) - y * Math.Sin(angle_RAD);
            //    histo_Y[i] = x * Math.Sin(angle_RAD) + y * Math.Cos(angle_RAD);
            //}

            my_form.processing();

            for (int i = 0; i < N; i++)
            {

                for (int k = 0; k < N; k++)
                {
                    double x = pointS_Hist[i,k].x;
                    double y = pointS_Hist[i,k].y;

                    double y_for_x = pointS_Hist[i, k].y;
                    pointS_Hist[i, k].x = x * Math.Cos(angle_RAD) - y * Math.Sin(angle_RAD);
                    pointS_Hist[i, k].y = x * Math.Sin(angle_RAD) + y * Math.Cos(angle_RAD);
                }
            }

            //for (int i = 0; i < pointS_Hist.Length; i++)
            //{
            //    for (int k = 0; k < pointS_Hist.Length; k++)
            //    {
            //        double x = pointS_Hist[i].x;
            //        double y = pointS_Hist[k].y;

            //        pointS_Hist[i].x = x * Math.Cos(angle_RAD) - y * Math.Sin(angle_RAD);
            //        pointS_Hist[k].y = x * Math.Sin(angle_RAD) + y * Math.Cos(angle_RAD);
            //    }
            //}
            my_form.hist_2d_comp(pointS_Hist, N);
        }







        static double[] mid_comp;
        static List<double>[] sheet_centered;
        public static void to_independent_nD(Form1 my_form)
        {
            mid_comp = new double[my_form.dimension];
            sheet_centered = new List<double>[my_form.dimension];
            for (int i = 0; i < my_form.dimension; i++)
            {
                sheet_centered[i] = new List<double>();
                mid_comp[i] = my_form.midarif_list[i]; 
                for (int k = 0; k < my_form.sheet[0].Count; k++)
                {
                    double not_ref = my_form.sheet[i][k] - my_form.midarif_list[i];
                    //double not_ref = my_form.sheet[i][k];
                    sheet_centered[i].Add(not_ref);
                }
            }
            ;
            List<double>[] trns = my_form.trans(sheet_centered, my_form.dimension, my_form.sheet[0].Count);
            double[][] new_DC_matrix = my_form.mult(sheet_centered, trns, my_form.dimension, my_form.sheet[0].Count, my_form.sheet[0].Count, my_form.dimension);
            double[,] matrix_arr = new double[my_form.dimension, my_form.dimension];
            for (int i = 0; i < my_form.dimension; i++)
            {
                for (int k = 0; k < my_form.dimension; k++)
                {
                    new_DC_matrix[i][k] /= my_form.sheet[0].Count;
                    matrix_arr[i, k] = new_DC_matrix[i][k];
                }
            }
            ;
            //Run_LU(new_DC_matrix, my_form.dimension);
            //Yakobi(my_form.dimension, new_DC_matrix

            Matrix<double> matrix = Matrix<double>.Build.DenseOfArray(matrix_arr);
            eigenVectors = eigenV(matrix);
            Evd<double> eigen = matrix.Evd();

            // Отримуємо власні вектори
            MathNet.Numerics.LinearAlgebra.Vector<Complex> eigenValues = eigen.EigenValues;// arr[]

            for (int i = 0; i < my_form.dimension; i++) //до n
            {
                for (int k = 0; k < my_form.sheet[0].Count; k++) // дo N
                {
                    double sum = 0;
                    for (int a = 0; a < my_form.dimension; a++)
                    {
                        sum += eigenVectors.Column(my_form.dimension - 1 - i)[a] * sheet_centered[a][k];
                       // sum += eigenVectors.Column(i)[a] * sheet_centered[a][k];
                    }
                    my_form.sheet[i][k] = sum;
                }
            }
            my_form.processing_multi();
            my_form.addlambda(eigenVectors, eigenValues); // зворотній порядок
                                                          //my_form.addlambda(roots, arrplus); // зворотній порядок
                                                          // my_form.addlambda(roots, lamDAs); // зворотній порядок


        }
        static Matrix<double> eigenVectors;
        static double[][] arrplus;
        public static void to_orig_nD(Form1 my_form)
        {
            char[] sep = new char[] { ' ' };
            string[] choosen = Prompt.ShowDialog("Головні компоненти", "").Split(sep, StringSplitOptions.RemoveEmptyEntries);
            int[] comp = new int[choosen.Length];
            for (int i = 0; i < choosen.Length; i++)
            {
                comp[i] = int.Parse(choosen[i]) - 1;
            }
            //if (choosen.Length > 2)
            //{
            double[,] temp_sheet = new double[my_form.dimension, my_form.sheet[0].Count];
            for (int i = 0; i < choosen.Length; i++)
            {
                for (int k = 0; k < my_form.sheet[0].Count; k++)
                {
                    double sum = 0;
                    for (int a = 0; a < my_form.dimension; a++)
                    {
                        sum += eigenVectors.Column(my_form.dimension - 1 - a)[comp[i]] * my_form.sheet[a][k];
                        //sum += eigenVectors.Column(a)[i] * my_form.sheet[a][k];
                    }
                    //my_form.sheet[comp[i]][k] = sum;
                    temp_sheet[comp[i], k] = sum;
                }

            }
            for (int i = 0; i < choosen.Length; i++)
            {
                for (int k = 0; k < my_form.sheet[0].Count; k++)
                {

                    my_form.sheet[comp[i]][k] = temp_sheet[comp[i], k];
                }

            }
            if(choosen.Length == my_form.dimension)
            {

                my_form.processing_multi();
            }
            else if (choosen.Length > 2)
            {

                //for (int i = 0; i < my_form.dimension; i++)
                //{
                //    for (int k = 0; k < my_form.sheet[0].Count; k++)
                //    {
                //        my_form.sheet[i][k] = my_form.sheet[i][k] + mid_comp[i];
                //        //double not_ref = my_form.sheet[i][k];
                //        //sheet_centered[i].Add(not_ref);
                //    }
                //}
                my_form.processing_multi();

            }
            else if (choosen.Length == 2)
            {
                my_form.vector1 = comp[0];
                my_form.vector2 = comp[1];
                my_form.all_2n = new List<double[]>();
                for (int i = 0; i < my_form.sheet[0].Count; i++)
                {
                    double[] temporary = new double[2];
                    temporary[0] = temp_sheet[my_form.vector1, i];
                    temporary[1] = temp_sheet[my_form.vector2, i];
                    my_form.all_2n.Add(temporary);
                    if (i == 0)
                    {
                        my_form.max_n2 = new double[2];
                        my_form.min_n2 = new double[2];
                        my_form.max_n2[0] = temp_sheet[i, my_form.vector1];
                        my_form.max_n2[1] = temp_sheet[i, my_form.vector2];
                        my_form.min_n2[0] = temp_sheet[i, my_form.vector1];
                        my_form.min_n2[1] = temp_sheet[i, my_form.vector2];
                    }
                    else
                    {
                        if (my_form.max_n2[0] < my_form.all_2n[i][my_form.vector1])
                        {
                            my_form.max_n2[0] = my_form.all_2n[i][my_form.vector1];
                        }
                        if (my_form.max_n2[1] < my_form.all_2n[i][my_form.vector2])
                        {
                            my_form.max_n2[1] = my_form.all_2n[i][my_form.vector2];
                        }
                        if (my_form.min_n2[0] > my_form.all_2n[i][my_form.vector1])
                        {
                            my_form.min_n2[0] = my_form.all_2n[i][my_form.vector1];
                        }
                        if (my_form.min_n2[1] > my_form.all_2n[i][my_form.vector2])
                        {
                            my_form.min_n2[1] = my_form.all_2n[i][my_form.vector2];
                        }
                    }

                }
                my_form.multi_space_sample = true;
                my_form.more_thn_2 = false;

                my_form.processing();
            }
            else if (choosen.Length == 1)
            {
                my_form.select_also(choosen[0]);
            }
            // }
        }



        public static MathNet.Numerics.LinearAlgebra.Matrix<double> eigenV(Matrix<double> matrix)
        {
            // Обчислюємо власні вектори
            Evd<double> eigen = matrix.Evd();

            // Отримуємо власні вектори
            MathNet.Numerics.LinearAlgebra.Vector<Complex> eigenValues = eigen.EigenValues;
            Matrix<double> eigenVectors = eigen.EigenVectors;

            return eigenVectors;
        }

        static void Yakobi(int n, double[][] arr1)
        {

            arr = new double[n][];// метод Якобі
            List<double[][]> arrays = new List<double[][]>();
            List<double[][]> arrays_U = new List<double[][]>();
            for (int h = 0; h < n; h++)
            {
                arr[h] = new double[n];
                for (int k = 0; k < n; k++)
                {
                    arr[h][k] = arr1[h][k];
                }
            }
            for (int u = 0; u < n; u++)
            {
                for (int k = 0; k < n; k++)
                {
                    if (arr[u][k] != arr[k][u])
                    {
                        MessageBox.Show("Матриця не є діагональною");
                    }
                }
            }
            double epsilon = 0.1;
            List<double[][]> multi_U = new List<double[][]>();
            double[][] E = new double[n][];
            for (int k = 0; k < n; k++)
            {
                E[k] = new double[n];
                for (int t = 0; t < n; t++)
                {
                    if (t == k)
                    {
                        E[k][t] = 1.0;
                    }
                    else
                    {
                        E[k][t] = 0.0;
                    }
                }
            }
            multi_U.Add(E);
            arrays.Add(arr);
            int i = 0;
            for (; ; i++)/////////////////////////////////////////////
            {
                List<double> for_max = new List<double>();
                double[][] A = new double[n][];
                int index = 0;
                for (int c = 0; c < n; c++)
                {
                    for (int k = 0; k < n; k++)
                    {

                        if (c != k)
                        {
                            for_max.Add(arrays[i][c][k]);
                            //for_max[index] = arrays[i][c][k];
                            index++;
                        }
                    }
                }
                double min = for_max.Min();
                double max = for_max.Max();
                double realmax;
                if (Math.Abs(Math.Round(min, 5)) > Math.Abs(Math.Round(max, 5)))
                {
                    realmax = min;
                }
                else
                {
                    realmax = max;
                }
                int a = 0;
                int b = 0;
                for (int t = 0; t < n; t++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        if (arrays[i][t][k] == realmax)
                        {
                            a = t;
                            b = k;
                            goto here;
                        }
                    }
                }
                here:
                ;
                double phi = 0.5 * Math.Atan(2 * realmax / (arrays[i][a][a] - arrays[i][b][b]));
                double[][] U = new double[n][];
                for (int k = 0; k < n; k++)
                {
                    U[k] = new double[n];
                    for (int r = 0; r < n; r++)
                    {
                        if (k == r)
                        {
                            U[k][r] = 1.0;
                        }
                        else
                        {
                            U[k][r] = 0;
                        }
                    }
                }
                U[a][a] = Math.Cos(phi);
                U[a][b] = Math.Sin(phi);
                U[b][a] = (-1) * Math.Sin(phi);
                U[b][b] = U[a][a];
                arrays_U.Add(U);
                double[][] multiplication = new double[n][];
                double[][] multiplication1 = new double[n][];
                double[][] trans = new double[n][];
                for (int f = 0; f < n; f++)
                {
                    multiplication[f] = new double[n];
                    multiplication1[f] = new double[n];
                    trans[f] = new double[n];
                    for (int o = 0; o < n; o++)
                    {
                        trans[f][o] = U[o][f];
                    }
                }
                for (int sho = 0; sho < n; sho++)
                {
                    for (int nicho = 0; nicho < n; nicho++)
                    {
                        for (int u = 0; u < n; u++)
                        {
                            multiplication[sho][nicho] = multiplication[sho][nicho] + U[sho][u] * arrays[i][u][nicho];
                        }
                    }
                }
                for (int sho = 0; sho < n; sho++)
                {
                    for (int nicho = 0; nicho < n; nicho++)
                    {
                        for (int u = 0; u < n; u++)
                        {
                            multiplication1[sho][nicho] = multiplication1[sho][nicho] + multiplication[sho][u] * trans[u][nicho];
                        }
                    }
                }
                arrays.Add(multiplication1);
                double sum = 0;
                for (int k = 0; k < n; k++)
                {
                    for (int v = 0; v < n; v++)
                    {
                        if (k != v)
                        {
                            sum = sum + Math.Pow(multiplication1[k][v], 2);
                        }
                    }
                }

                double[][] temp = new double[n][];
                for (int k = 0; k < n; k++)
                {
                    temp[k] = new double[n];
                }
                for (int s = 0; s < n; s++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        for (int r = 0; r < n; r++)
                        {
                            temp[s][k] = temp[s][k] + multi_U[i][s][r] * arrays_U[i][r][k];
                        }
                    }
                }
                multi_U.Add(temp);

                //MessageBox.Show(sum.ToString());
                if (sum <= epsilon)
                {
                    i++;
                    break;
                }


            }

            //richTextBox1.Text = "Власні числа:\n";
            lamDAs = new double[n];
            for (int o = 0; o < n; o++)
            {
                for (int k = 0; k < n; k++)
                {
                    if (o == k)
                    {
                        lamDAs[k] = arrays[i][o][k];
                        //richTextBox1.Text = richTextBox1.Text + Math.Round(arrays[i][o][k], 5).ToString() + "\n";
                    }
                }
            }
            double[][][] resulting = new double[i][][];
            for (int t = 0; t < i; t++)
            {
                resulting[t] = new double[n][];
                for (int a = 0; a < n; a++)
                {
                    resulting[t][a] = new double[n];
                }
            }

            // richTextBox1.Text = richTextBox1.Text + "Власні вектори:\n";
            roots = new double[n][];
            for (int a = 0; a < n; a++)
            {
                roots[a] = new double[n];
                //richTextBox1.Text = richTextBox1.Text + "(";
                for (int b = 0; b < n; b++)
                {
                    roots[a][b] = multi_U[i][b][a];
                    //if (b == (n - 1))
                    //{
                    //    richTextBox1.Text = richTextBox1.Text + Math.Round(multi_U[i][b][a], 5).ToString();
                    //}
                    //else
                    //{
                    //    richTextBox1.Text = richTextBox1.Text + Math.Round(multi_U[i][b][a], 5).ToString() + "; ";
                    //}
                }
                //richTextBox1.Text = richTextBox1.Text + ")\n";
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //double[] X_0 = new double[n];//метод степеневий
            //X_0[0] = 1.0;
            //for (int k = 1; k < n; k++)
            //{
            //    X_0[k] = 1.0;
            //    }
            //    List<double[]> X = new List<double[]>();
            //    List<double> lambda_max = new List<double>();
            //    X.Add(X_0);
            //    richTextBox1.Text = richTextBox1.Text + "Максимальне власне значення:\n";
            //    for (int k = 0; k < n; k++)
            //    {
            //        double[] tem = new double[n];


            //        for (int s = 0; s < n; s++)
            //        {
            //            for (int a = 0; a < n; a++)
            //            {
            //                tem[s] = tem[s] + arr[s][a] * X[k][a];
            //            }
            //        }
            //        X.Add(tem);
            //        double norm1 = 0;
            //        double norm2 = 0;
            //        for (int g = 0; g < n; g++)
            //        {
            //            norm1 = norm1 + Math.Pow(X[k + 1][g], 2);
            //            norm2 = norm2 + Math.Pow(X[k][g], 2);
            //        }
            //        double lambdaa = Math.Sqrt(norm1 / norm2);
            //        lambda_max.Add(lambdaa);
            //    }
            //    richTextBox1.Text = richTextBox1.Text + Math.Round(lambda_max[n - 1], 5).ToString() + "\n";
            //    double norm;
            //    List<double> l = new List<double>();
            //    for (int k = 0; k < n; k++)
            //    {
            //        l.Add(Math.Abs(X[n][k]));
            //    }
            //    norm = l.Max();
            //    double[] max_vector = new double[n];
            //    for (int k = 0; k < n; k++)
            //    {
            //        max_vector[k] = X[n][k] / norm;
            //    }
            //    richTextBox1.Text = richTextBox1.Text + "Максимальний вектор:\n";
            //    for (int k = 0; k < n; k++)
            //    {
            //        richTextBox1.Text = richTextBox1.Text + $"{max_vector[k]}\n";
            //    }


            //    richTextBox1.Text = richTextBox1.Text + "Мінімальне власне значення:\n";
            //    List<double> lambda_min = new List<double>();
            //    X = new List<double[]>();
            //    X_0 = new double[n];
            //    X_0[0] = 1.0;
            //    for (int k = 1; k < n; k++)
            //    {
            //        X_0[k] = 0.0;
            //    }

            //    X.Add(X_0);
            //    double[][] inverteed = inverted();
            //    for (int k = 0; k < n; k++)
            //    {
            //        double[] tem = new double[n];

            //        for (int s = 0; s < n; s++)
            //        {
            //            for (int a = 0; a < n; a++)
            //            {
            //                tem[s] = tem[s] + inverteed[s][a] * X[k][a];
            //            }
            //        }
            //        X.Add(tem);
            //        double norm1 = 0;
            //        double norm2 = 0;
            //        for (int g = 0; g < n; g++)
            //        {
            //            norm1 = norm1 + Math.Pow(X[k + 1][g], 2);
            //            norm2 = norm2 + Math.Pow(X[k][g], 2);
            //        }
            //        double lambdaa = Math.Sqrt(norm2 / norm1);
            //        lambda_min.Add(lambdaa);
            //    }
            //    richTextBox1.Text = richTextBox1.Text + Math.Round(lambda_min[n - 1], 5).ToString() + "\n";
            //    norm = 0;
            //    l = new List<double>();
            //    for (int k = 0; k < n; k++)
            //    {
            //        l.Add(Math.Abs(X[n][k]));
            //    }
            //    norm = l.Max();
            //    max_vector = new double[n];
            //    for (int k = 0; k < n; k++)
            //    {
            //        max_vector[k] = X[n][k] / norm;
            //    }
            //    richTextBox1.Text = richTextBox1.Text + "Мінімальний вектор:\n";
            //    for (int k = 0; k < n; k++)
            //    {
            //        richTextBox1.Text = richTextBox1.Text + $"{max_vector[k]}\n";
            //    }
            //    richTextBox1.Text = richTextBox1.Text + $"Кількість ітерацій {i + 1}\n";

        }



        //LU below
        int count = 0;
        static double[][] arr;
        static double[][] L;
        static double[][] U;
        static double[][] roots;
        static double[] lamDAs;
        //static void LU(double[][] array, int n)
        //{
        //    if (array[0][0] == 0)
        //    {
        //        MessageBox.Show("A11 == 0");
        //    }
        //    else
        //    {
        //        //make();
        //        U = new double[n][];
        //        L = new double[n][];
        //        for (int i = 0; i < n; i++)
        //        {
        //            U[i] = new double[n];
        //            L[i] = new double[n];
        //        }
        //        for (int i = 0; i < n; i++)
        //        {
        //            double sum = 0;
        //            for (int k = 0; k < n; k++)
        //            {
        //                sum = 0;
        //                if (i >= k)
        //                {
        //                    for (int l = 0; l <= k - 1; l++)
        //                    {
        //                        sum = sum + L[i][l] * U[l][k];
        //                    }
        //                    L[i][k] = array[i][k] - sum;
        //                    //arr3[i][k].Text = L[i][k].ToString();
        //                }
        //                if (i == 1 && k == 2)
        //                {
        //                    ;
        //                }
        //                if (i < k)
        //                {
        //                    sum = 0;
        //                    for (int l = 0; l <= i - 1; l++)
        //                    {
        //                        sum = sum + L[i][l] * U[l][k];
        //                    }
        //                    U[i][k] = (array[i][k] - sum) / L[i][i];
        //                    //arr5[i][k].Text = U[i][k].ToString();
        //                }
        //            }
        //            U[i][i] = 1;
        //            //arr5[i][i].Text = U[i][i].ToString();
        //        }
        //    }
        //}
        //private static void Run_LU(double[][]DC, int n)
        //{
        //    roots = new double[n][];
        //    arr = new double[n][];
        //    arrplus = new double[n][];
        //    for (int i = 0; i < n; i++)
        //    {
        //        roots[i] = new double[n];
        //        arr[i] = new double[n];
        //        arrplus[i] = new double[n];
        //        for (int k = 0; k < n; k++)
        //        {
        //            arr[i][k] = DC[i][k];
        //        }
        //    }
        //    for (; ; )
        //    {
        //        LU(arr, n);
        //        for (int i = 0; i < n; i++)
        //        {
        //            for (int k = 0; k < n; k++)
        //            {
        //                arrplus[i][k] = 0;
        //            }
        //        }
        //        for (int i = 0; i < n; i++)
        //        {
        //            for (int k = 0; k < n; k++)
        //            {
        //                for (int a = 0; a < n; a++)
        //                {
        //                    arrplus[i][k] = arrplus[i][k] + U[i][a] * L[a][k];
        //                }
        //            }
        //        }
        //        ;
        //        int count = 0;
        //        for (int i = 0; i < n; i++)
        //        {
        //            if (Math.Abs(arrplus[i][i] - arr[i][i]) <= 0.0001)
        //            {
        //                count++;
        //            }
        //        }
        //        if (count == n)
        //        {
        //            goto here;
        //        }
        //        for (int i = 0; i < n; i++)
        //        {
        //            for (int k = 0; k < n; k++)
        //            {
        //                arr[i][k] = arrplus[i][k];
        //            }
        //        }
        //    }
        //    here:
        //    ;
        //    //richTextBox1.Text = "";
        //    double[] lambda = new double[n];
        //    for (int i = 0; i < n; i++)
        //    {
        //      //  richTextBox1.Text = richTextBox1.Text + arrplus[i][i].ToString() + "\n";
        //        lambda[i] = arrplus[i][i];
        //    }
        //    for (int i = 0; i < n; i++)
        //    {
        //        for (int k = 0; k < n; k++)
        //        {
        //            arr[i][k] = DC[i][k];
        //        }
        //    }
        //    for (int i = 0; i < n; i++)
        //    {
        //        double[][] for_solve = new double[n][];
        //        for (int k = 0; k < n; k++)
        //        {
        //            for_solve[k] = new double[n];
        //            for (int a = 0; a < n; a++)
        //            {
        //                for_solve[k][a] = arr[k][a];
        //            }
        //        }
        //        for (int k = 0; k < n; k++)
        //        {
        //            for_solve[k][k] = for_solve[k][k] - lambda[i];
        //        }
        //        ;
        //        for (int k = 0; k < n; k++)
        //        {
        //            for (int a = k + 1; a < n; a++)
        //            {
        //                double t = for_solve[a][k];
        //                for (int b = 0; b < n; b++)
        //                {
        //                    for_solve[a][b] = for_solve[a][b] - for_solve[k][b] * t / for_solve[k][k];
        //                }
        //            }
        //        }
        //        //richTextBox1.Text = richTextBox1.Text + "(";
        //        double[] b_arr = new double[n];
        //        for (int k = 0; k < n; k++)
        //        {
        //            b_arr[k] = (-1) * for_solve[k][n - 1];
        //        }

        //        roots[i][n - 1] = 1;
        //        for (int k = n - 2; k >= 0; k--)
        //        {
        //            double sum = 0;
        //            for (int a = k; a <= n - 2; a++)
        //            {
        //                sum = sum + for_solve[k][a] * roots[i][a];
        //            }
        //            roots[i][k] = (b_arr[k] - sum) / for_solve[k][k];
        //        }
        //        //for (int k = 0; k < n; k++)
        //        //{
        //        //    if (k == n - 1)
        //        //    {
        //        //        richTextBox1.Text = richTextBox1.Text + roots[k].ToString() + ")";
        //        //    }
        //        //    else
        //        //    {
        //        //        richTextBox1.Text = richTextBox1.Text + roots[k].ToString() + "; ";
        //        //    }
        //        //}
        //        ;
        //    }
        //    ;
        //}
    }
}

