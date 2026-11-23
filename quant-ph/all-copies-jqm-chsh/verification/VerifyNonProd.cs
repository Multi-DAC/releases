// Exact global verifier over all 55,720,396,530 spanned 4-wing subsets for the
// non-product two-copy form K_np = D_t (x) D_t + D_Delta (x) D_H (Proposition 11).
// Step 1: python -I -B -X utf8 prep_nonprod_exact.py   (writes nonprod_data.bin)
// Step 2: C:\Windows\Microsoft.NET\Framework64\v4.0.30319\csc.exe /optimize+ /platform:x64 /r:System.Numerics.dll VerifyNonProd.cs
// Step 3: VerifyNonProd.exe [nonprod_data.bin]
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;

class VerifyNonProd
{
    static long[] Zt = new long[256];
    static long[] ZDelta = new long[256];
    static long[] ZHnull = new long[256];
    static long qt_val, qd_val, qh_val;
    static BigInteger qt2, scaleDelta;

    static double[] B3 = new double[64 * 27];
    static double[] A00 = new double[27 * 27];
    static double[] A33 = new double[27 * 27];
    static double[] A03 = new double[27 * 27];
    static double[] B11 = new double[27 * 27];
    static double[] B22 = new double[27 * 27];
    static double[] B12 = new double[27 * 27];
    static double[] Xmap = new double[18 * 27];
    static double[] Ymap = new double[18 * 27];

    static ushort[] S2;
    static ulong[] S3;
    static int[] S3_s1, S3_s2;

    static double[] Z3, A03Z3, B12Z3, X3, Y3;
    static double[] a00v, a33v, b11v, b22v;

    const double PRUNE_THRESH = 1e-5;

    static void LoadData(string path)
    {
        using (var br = new BinaryReader(File.OpenRead(path)))
        {
            for (int i = 0; i < 256; i++) Zt[i] = br.ReadInt64();
            for (int i = 0; i < 256; i++) ZDelta[i] = br.ReadInt64();
            for (int i = 0; i < 256; i++) ZHnull[i] = br.ReadInt64();
            qt_val = br.ReadInt64();
            qd_val = br.ReadInt64();
            qh_val = br.ReadInt64();
            for (int i = 0; i < 64 * 27; i++) B3[i] = br.ReadDouble();
            for (int i = 0; i < 27 * 27; i++) A00[i] = br.ReadDouble();
            for (int i = 0; i < 27 * 27; i++) A33[i] = br.ReadDouble();
            for (int i = 0; i < 27 * 27; i++) A03[i] = br.ReadDouble();
            for (int i = 0; i < 27 * 27; i++) B11[i] = br.ReadDouble();
            for (int i = 0; i < 27 * 27; i++) B22[i] = br.ReadDouble();
            for (int i = 0; i < 27 * 27; i++) B12[i] = br.ReadDouble();
            for (int i = 0; i < 18 * 27; i++) Xmap[i] = br.ReadDouble();
            for (int i = 0; i < 18 * 27; i++) Ymap[i] = br.ReadDouble();
        }
        qt2 = (BigInteger)qt_val * (BigInteger)qt_val;
        BigInteger qdh = (BigInteger)qd_val * (BigInteger)qh_val;
        BigInteger rem;
        scaleDelta = BigInteger.DivRem(qt2, qdh, out rem);
        if (!rem.IsZero) throw new Exception("qt^2 not divisible by q_Delta * q_Hnull!");
    }

    static void BuildSpannedSets()
    {
        var s1 = new HashSet<int>();
        for (int m = 0; m < 16; m++)
            if (((m >> 0) & 1) + ((m >> 3) & 1) == ((m >> 1) & 1) + ((m >> 2) & 1))
                s1.Add(m);

        var s2List = new List<ushort>();
        for (int m = 0; m < 65536; m++)
        {
            bool ok = true;
            for (int wa = 0; wa < 4; wa++)
                if (!s1.Contains((m >> (4 * wa)) & 15)) { ok = false; break; }
            if (!ok) continue;
            for (int wb = 0; wb < 4; wb++)
            {
                int col = ((m >> wb) & 1) | (((m >> (4 + wb)) & 1) << 1)
                        | (((m >> (8 + wb)) & 1) << 2) | (((m >> (12 + wb)) & 1) << 3);
                if (!s1.Contains(col)) { ok = false; break; }
            }
            if (ok) s2List.Add((ushort)m);
        }
        S2 = s2List.ToArray();
        if (S2.Length != 90) throw new Exception("Expected |S2|=90");

        var pairs2 = new Dictionary<uint, List<uint>>();
        for (int i = 0; i < 90; i++)
        {
            uint a = S2[i];
            for (int j = 0; j < 90; j++)
            {
                uint d = S2[j];
                uint key = (a ^ d) | ((a & d) << 16);
                List<uint> lst;
                if (!pairs2.TryGetValue(key, out lst)) { lst = new List<uint>(); pairs2[key] = lst; }
                lst.Add(a | (d << 16));
            }
        }
        var s3List = new List<ulong>(51678);
        foreach (var lst in pairs2.Values)
        {
            int n = lst.Count;
            for (int i = 0; i < n; i++)
            {
                ulong a = lst[i] & 0xFFFFUL;
                ulong d = (lst[i] >> 16) & 0xFFFFUL;
                for (int j = 0; j < n; j++)
                {
                    ulong b = lst[j] & 0xFFFFUL;
                    ulong c = (lst[j] >> 16) & 0xFFFFUL;
                    s3List.Add(a | (b << 16) | (c << 32) | (d << 48));
                }
            }
        }
        s3List.Sort();
        S3 = s3List.ToArray();
        if (S3.Length != 51678) throw new Exception("Expected |S3|=51678");
        S3_s1 = new int[S3.Length];
        S3_s2 = new int[S3.Length];
        for (int i = 0; i < S3.Length; i++)
        {
            S3_s1[i] = (int)((S3[i] >> 16) & 0xFFFFUL);
            S3_s2[i] = (int)((S3[i] >> 32) & 0xFFFFUL);
        }
    }

    static void PrecomputeS3()
    {
        int N = S3.Length;
        Z3 = new double[N * 27];
        A03Z3 = new double[N * 27];
        B12Z3 = new double[N * 27];
        X3 = new double[N * 18];
        Y3 = new double[N * 18];
        a00v = new double[N];
        a33v = new double[N];
        b11v = new double[N];
        b22v = new double[N];

        Parallel.For(0, N, idx =>
        {
            ulong m = S3[idx];
            int zOff = idx * 27;
            ulong tmp = m;
            while (tmp != 0UL)
            {
                ulong lsb = tmp & (0UL - tmp);
                int bit = TrailingZeroCount(lsb);
                int bOff = bit * 27;
                for (int k = 0; k < 27; k++) Z3[zOff + k] += B3[bOff + k];
                tmp ^= lsb;
            }
            double sA00 = 0, sA33 = 0, sB11 = 0, sB22 = 0;
            for (int r = 0; r < 27; r++)
            {
                double zr = Z3[zOff + r];
                int rOff = r * 27;
                double dA00 = 0, dA33 = 0, dA03 = 0, dB11 = 0, dB22 = 0, dB12 = 0;
                for (int c = 0; c < 27; c++)
                {
                    double zc = Z3[zOff + c];
                    dA00 += A00[rOff + c] * zc;
                    dA33 += A33[rOff + c] * zc;
                    dA03 += A03[rOff + c] * zc;
                    dB11 += B11[rOff + c] * zc;
                    dB22 += B22[rOff + c] * zc;
                    dB12 += B12[rOff + c] * zc;
                }
                sA00 += zr * dA00;
                sA33 += zr * dA33;
                sB11 += zr * dB11;
                sB22 += zr * dB22;
                A03Z3[zOff + r] = dA03;
                B12Z3[zOff + r] = dB12;
            }
            a00v[idx] = sA00;
            a33v[idx] = sA33;
            b11v[idx] = sB11;
            b22v[idx] = sB22;

            int xOff = idx * 18;
            for (int r = 0; r < 18; r++)
            {
                int rOff = r * 27;
                double dx = 0, dy = 0;
                for (int c = 0; c < 27; c++)
                {
                    double zc = Z3[zOff + c];
                    dx += Xmap[rOff + c] * zc;
                    dy += Ymap[rOff + c] * zc;
                }
                X3[xOff + r] = dx;
                Y3[xOff + r] = dy;
            }
        });
    }

    static int TrailingZeroCount(ulong x)
    {
        int c = 0;
        if ((x & 0xFFFFFFFFUL) == 0) { x >>= 32; c += 32; }
        if ((x & 0xFFFFUL) == 0) { x >>= 16; c += 16; }
        if ((x & 0xFFUL) == 0) { x >>= 8; c += 8; }
        if ((x & 0xFUL) == 0) { x >>= 4; c += 4; }
        if ((x & 0x3UL) == 0) { x >>= 2; c += 2; }
        if ((x & 0x1UL) == 0) { c += 1; }
        return c;
    }

    static BigInteger EvalTerm(long[] M1, long[] M2, int[] F)
    {
        BigInteger[] MF = new BigInteger[256];
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 16; j++)
            {
                BigInteger s = BigInteger.Zero;
                for (int k = 0; k < 16; k++)
                    if (F[16 * k + j] != 0) s += M1[16 * i + k];
                MF[16 * i + j] = s;
            }
        }
        BigInteger total = BigInteger.Zero;
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 16; j++)
            {
                if (F[16 * i + j] != 0)
                {
                    BigInteger row = BigInteger.Zero;
                    for (int l = 0; l < 16; l++)
                        row += MF[16 * i + l] * M2[16 * l + j];
                    total += row;
                }
            }
        }
        return total;
    }

    static BigInteger EvalExactBigInt(ulong g0, ulong g1, ulong g2, ulong g3)
    {
        ulong[] gs = new ulong[] { g0, g1, g2, g3 };
        int[] F = new int[256];
        for (int w4 = 0; w4 < 4; w4++)
        {
            ulong g = gs[w4];
            for (int w3 = 0; w3 < 4; w3++)
            {
                int slice = (int)((g >> (16 * w3)) & 0xFFFFUL);
                int c2 = 4 * w3 + w4;
                for (int c1 = 0; c1 < 16; c1++)
                    F[16 * c1 + c2] = (slice >> c1) & 1;
            }
        }
        BigInteger t1 = EvalTerm(Zt, Zt, F);
        BigInteger t2 = EvalTerm(ZDelta, ZHnull, F);
        return t1 + scaleDelta * t2;
    }

    static bool IsSpanned4(ulong g0, ulong g1, ulong g2, ulong g3)
    {
        return ((g0 ^ g3) == (g1 ^ g2)) && ((g0 & g3) == (g1 & g2));
    }

    static double EvalFormula(int i0, int i1, int i2, int i3)
    {
        int o0 = i0 * 27, o3 = i3 * 27, o1 = i1 * 27, o2 = i2 * 27;
        double dotA = 0, dotB = 0;
        for (int k = 0; k < 27; k++)
        {
            dotA += Z3[o0 + k] * A03Z3[o3 + k];
            dotB += Z3[o1 + k] * B12Z3[o2 + k];
        }
        double alpha = a00v[i0] + a33v[i3] + dotA;
        double beta = b11v[i1] + b22v[i2] + dotB;
        int x0 = i0 * 18, x3 = i3 * 18, y1 = i1 * 18, y2 = i2 * 18;
        double cross = 0;
        for (int k = 0; k < 18; k++)
            cross += (X3[x0 + k] - X3[x3 + k]) * (Y3[y1 + k] - Y3[y2 + k]);
        return alpha + beta + cross;
    }

    class S0Task
    {
        public uint S0Key;
        public List<int> PairsA, PairsD;
        public long EstPairs;
    }

    class MinRecord
    {
        public BigInteger ExactVal;
        public ulong g0, g1, g2, g3;
    }

    static void Main(string[] args)
    {
        var sw = Stopwatch.StartNew();
        string dataPath = args.Length > 0 ? args[0] : Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "nonprod_data.bin");
        LoadData(dataPath);
        BuildSpannedSets();
        PrecomputeS3();
        Console.WriteLine("Precomputed S3 (51678) in {0:F2}s", sw.Elapsed.TotalSeconds);

        double qt2d = (double)qt2;
        var rng = new Random(54321);
        double maxSelfErr = 0.0;
        var s3Index = new Dictionary<ulong, int>(51678);
        for (int i = 0; i < S3.Length; i++) s3Index[S3[i]] = i;

        int tested = 0;
        while (tested < 500)
        {
            int i0 = rng.Next(S3.Length);
            int i1 = rng.Next(S3.Length);
            ulong g0 = S3[i0], g1 = S3[i1];
            ulong h1 = g1 ^ g0;
            int start = rng.Next(S3.Length);
            int i2 = -1;
            for (int t = 0; t < 200; t++)
            {
                int cand = (start + t) % S3.Length;
                if (((S3[cand] ^ g0) & h1) == 0UL) { i2 = cand; break; }
            }
            if (i2 < 0) continue;
            ulong g2 = S3[i2];
            ulong g3 = g1 ^ g2 ^ g0;
            int i3 = s3Index[g3];
            double fVal = EvalFormula(i0, i1, i2, i3);
            BigInteger exVal = EvalExactBigInt(g0, g1, g2, g3);
            double exD = (double)exVal / qt2d;
            double err = Math.Abs(fVal - exD);
            if (err > maxSelfErr) maxSelfErr = err;
            tested++;
        }
        Console.WriteLine("Self-test (500 random spanned quadruples vs exact BigInteger): max_err = {0:E3}", maxSelfErr);
        if (maxSelfErr > 1e-12) throw new Exception("Self-test error too large!");

        var s2IdxMap = new Dictionary<ushort, int>();
        for (int i = 0; i < 90; i++) s2IdxMap[S2[i]] = i;
        var byA0List = new List<int>[90];
        for (int i = 0; i < 90; i++) byA0List[i] = new List<int>();
        for (int i = 0; i < S3.Length; i++)
            byA0List[s2IdxMap[(ushort)(S3[i] & 0xFFFFUL)]].Add(i);
        var byA0 = new int[90][];
        for (int i = 0; i < 90; i++) byA0[i] = byA0List[i].ToArray();

        var s0Dict = new Dictionary<uint, S0Task>();
        for (int a = 0; a < 90; a++)
        {
            uint av = S2[a];
            for (int d = a; d < 90; d++)
            {
                uint dv = S2[d];
                uint key = (av ^ dv) | ((av & dv) << 16);
                S0Task task;
                if (!s0Dict.TryGetValue(key, out task))
                {
                    task = new S0Task { S0Key = key, PairsA = new List<int>(), PairsD = new List<int>(), EstPairs = 0 };
                    s0Dict[key] = task;
                }
                task.PairsA.Add(a);
                task.PairsD.Add(d);
                long na = byA0[a].Length, nd = byA0[d].Length;
                task.EstPairs += (a == d) ? (na * (na + 1) / 2) : (na * nd);
            }
        }
        var tasks = new List<S0Task>(s0Dict.Values);
        tasks.Sort((x, y) => y.EstPairs.CompareTo(x.EstPairs));

        long totalQuadruples = 0, exactChecked = 0, negFound = 0;
        object minLock = new object();
        BigInteger globalMinExact = BigInteger.Pow(10, 50);
        var globalMinimizers = new List<MinRecord>();
        int doneCount = 0;
        var swRun = Stopwatch.StartNew();
        var opts = new ParallelOptions { MaxDegreeOfParallelism = 16 };

        Parallel.For(0, tasks.Count, opts, tIdx =>
        {
            var task = tasks[tIdx];
            int cap = (int)task.EstPairs;
            ulong[] keys = new ulong[cap];
            uint[] pairs = new uint[cap];
            int cnt = 0;

            for (int p = 0; p < task.PairsA.Count; p++)
            {
                int aIdx = task.PairsA[p], dIdx = task.PairsD[p];
                int[] listA = byA0[aIdx], listD = byA0[dIdx];
                int nA = listA.Length, nD = listD.Length;
                if (aIdx == dIdx)
                {
                    for (int ia = 0; ia < nA; ia++)
                    {
                        int g0 = listA[ia], a1 = S3_s1[g0], a2 = S3_s2[g0];
                        for (int id = ia; id < nA; id++)
                        {
                            int g3 = listA[id], d1 = S3_s1[g3], d2 = S3_s2[g3];
                            keys[cnt] = (ulong)(uint)(a1 ^ d1) | ((ulong)(uint)(a1 & d1) << 16)
                                      | ((ulong)(uint)(a2 ^ d2) << 32) | ((ulong)(uint)(a2 & d2) << 48);
                            pairs[cnt] = (uint)g0 | ((uint)g3 << 16);
                            cnt++;
                        }
                    }
                }
                else
                {
                    for (int ia = 0; ia < nA; ia++)
                    {
                        int g0 = listA[ia], a1 = S3_s1[g0], a2 = S3_s2[g0];
                        for (int id = 0; id < nD; id++)
                        {
                            int g3 = listD[id], d1 = S3_s1[g3], d2 = S3_s2[g3];
                            keys[cnt] = (ulong)(uint)(a1 ^ d1) | ((ulong)(uint)(a1 & d1) << 16)
                                      | ((ulong)(uint)(a2 ^ d2) << 32) | ((ulong)(uint)(a2 & d2) << 48);
                            pairs[cnt] = (g0 < g3) ? ((uint)g0 | ((uint)g3 << 16)) : ((uint)g3 | ((uint)g0 << 16));
                            cnt++;
                        }
                    }
                }
            }

            Array.Sort(keys, pairs, 0, cnt);

            long locQuads = 0, locExact = 0, locNeg = 0;
            BigInteger locMinExact = BigInteger.Pow(10, 50);
            var locMinimizers = new List<MinRecord>();

            int maxM = 1024;
            int[] ordG0 = new int[maxM], ordG3 = new int[maxM];
            double[] alphaArr = new double[maxM], betaArr = new double[maxM], normArr = new double[maxM];
            double[] Xbuf = new double[maxM * 18], Ybuf = new double[maxM * 18];

            int pos = 0;
            while (pos < cnt)
            {
                ulong curKey = keys[pos];
                int end = pos + 1;
                while (end < cnt && keys[end] == curKey) end++;

                int m = 0;
                for (int q = pos; q < end; q++)
                {
                    uint pk = pairs[q];
                    m += ((pk & 0xFFFFU) == (pk >> 16)) ? 1 : 2;
                }
                locQuads += (long)m * m;

                if (m > maxM)
                {
                    maxM = m * 2;
                    ordG0 = new int[maxM]; ordG3 = new int[maxM];
                    alphaArr = new double[maxM]; betaArr = new double[maxM]; normArr = new double[maxM];
                    Xbuf = new double[maxM * 18]; Ybuf = new double[maxM * 18];
                }

                double minA = double.PositiveInfinity, minB = double.PositiveInfinity, maxN = 0.0;
                int idx = 0;
                for (int q = pos; q < end; q++)
                {
                    uint pk = pairs[q];
                    int g0 = (int)(pk & 0xFFFFU), g3 = (int)(pk >> 16);
                    int o0 = g0 * 27, o3 = g3 * 27;
                    double dotA = 0.0, dotB = 0.0;
                    for (int k = 0; k < 27; k++)
                    {
                        double z0k = Z3[o0 + k];
                        dotA += z0k * A03Z3[o3 + k];
                        dotB += z0k * B12Z3[o3 + k];
                    }
                    double a03 = a00v[g0] + a33v[g3] + dotA;
                    double b03 = b11v[g0] + b22v[g3] + dotB;

                    if (g0 == g3)
                    {
                        ordG0[idx] = g0; ordG3[idx] = g3;
                        alphaArr[idx] = a03; betaArr[idx] = b03; normArr[idx] = 0.0;
                        int xOff = idx * 18;
                        for (int k = 0; k < 18; k++) { Xbuf[xOff + k] = 0.0; Ybuf[xOff + k] = 0.0; }
                        if (a03 < minA) minA = a03;
                        if (b03 < minB) minB = b03;
                        idx++;
                    }
                    else
                    {
                        double a30 = a00v[g3] + a33v[g0] + dotA;
                        double b30 = b11v[g3] + b22v[g0] + dotB;
                        int x0 = g0 * 18, x3 = g3 * 18;
                        int xOff1 = idx * 18, xOff2 = (idx + 1) * 18;
                        double n2 = 0.0;
                        for (int k = 0; k < 18; k++)
                        {
                            double dx = X3[x0 + k] - X3[x3 + k];
                            double dy = Y3[x0 + k] - Y3[x3 + k];
                            Xbuf[xOff1 + k] = dx; Ybuf[xOff1 + k] = dy;
                            Xbuf[xOff2 + k] = -dx; Ybuf[xOff2 + k] = -dy;
                            n2 += dx * dx;
                        }
                        double nrm = Math.Sqrt(n2);
                        ordG0[idx] = g0; ordG3[idx] = g3; alphaArr[idx] = a03; betaArr[idx] = b03; normArr[idx] = nrm;
                        ordG0[idx + 1] = g3; ordG3[idx + 1] = g0; alphaArr[idx + 1] = a30; betaArr[idx + 1] = b30; normArr[idx + 1] = nrm;
                        if (a03 < minA) minA = a03;
                        if (a30 < minA) minA = a30;
                        if (b03 < minB) minB = b03;
                        if (b30 < minB) minB = b30;
                        if (nrm > maxN) maxN = nrm;
                        idx += 2;
                    }
                }

                if (minA + minB - maxN * maxN > PRUNE_THRESH) { pos = end; continue; }

                for (int i = 0; i < m; i++)
                {
                    double ai = alphaArr[i], ni = normArr[i];
                    if (ai + minB - ni * maxN > PRUNE_THRESH) continue;
                    int xOff = i * 18, gi0 = ordG0[i], gi3 = ordG3[i];
                    for (int j = 0; j < m; j++)
                    {
                        double bj = betaArr[j], nj = normArr[j];
                        if (ai + bj - ni * nj > PRUNE_THRESH) continue;
                        int yOff = j * 18;
                        double cross = 0.0;
                        for (int k = 0; k < 18; k++) cross += Xbuf[xOff + k] * Ybuf[yOff + k];
                        if (ai + bj + cross < PRUNE_THRESH)
                        {
                            int gj1 = ordG0[j], gj2 = ordG3[j];
                            if (gi0 == 0 && gj1 == 0 && gj2 == 0 && gi3 == 0) continue;
                            locExact++;
                            BigInteger ex = EvalExactBigInt(S3[gi0], S3[gj1], S3[gj2], S3[gi3]);
                            if (ex < 0)
                            {
                                locNeg++;
                                if (locNeg <= 10)
                                    Console.WriteLine("  NEG: g0={0:X16} g1={1:X16} g2={2:X16} g3={3:X16} val={4:E6}",
                                        S3[gi0], S3[gj1], S3[gj2], S3[gi3], (double)ex / qt2d);
                            }
                            if (ex < locMinExact)
                            {
                                locMinExact = ex;
                                locMinimizers.Clear();
                                locMinimizers.Add(new MinRecord { ExactVal = ex, g0 = S3[gi0], g1 = S3[gj1], g2 = S3[gj2], g3 = S3[gi3] });
                            }
                            else if (ex == locMinExact && locMinimizers.Count < 32)
                            {
                                locMinimizers.Add(new MinRecord { ExactVal = ex, g0 = S3[gi0], g1 = S3[gj1], g2 = S3[gj2], g3 = S3[gi3] });
                            }
                        }
                    }
                }
                pos = end;
            }

            Interlocked.Add(ref totalQuadruples, locQuads);
            Interlocked.Add(ref exactChecked, locExact);
            Interlocked.Add(ref negFound, locNeg);
            lock (minLock)
            {
                if (locMinExact < globalMinExact)
                {
                    globalMinExact = locMinExact;
                    globalMinimizers.Clear();
                    globalMinimizers.AddRange(locMinimizers);
                }
                else if (locMinExact == globalMinExact)
                {
                    foreach (var mr in locMinimizers)
                        if (globalMinimizers.Count < 64) globalMinimizers.Add(mr);
                }
            }
            int dNow = Interlocked.Increment(ref doneCount);
            if (dNow % 500 == 0 || dNow == tasks.Count)
            {
                Console.WriteLine("  [{0,6:F1}s] {1}/{2} done; quads={3}, exactChecked={4}, neg={5}",
                    swRun.Elapsed.TotalSeconds, dNow, tasks.Count,
                    Interlocked.Read(ref totalQuadruples), Interlocked.Read(ref exactChecked), Interlocked.Read(ref negFound));
            }
        });

        Console.WriteLine("\n======== EXACT GLOBAL NON-PRODUCT 2-COPY AUDIT COMPLETE ========");
        Console.WriteLine("Total time: {0:F2}s", sw.Elapsed.TotalSeconds);
        Console.WriteLine("Total spanned 4-wing subsets (incl. empty): {0} (expected 55720396530)", totalQuadruples);
        Console.WriteLine("Exact BigInteger checked (< 1e-5): {0}, negative count: {1}", exactChecked, negFound);
        Console.WriteLine("GLOBAL EXACT MIN q_t^2 * K(f, f) = {0}", globalMinExact);
        Console.WriteLine("GLOBAL EXACT MIN K(f, f) float   = {0:E15}", (double)globalMinExact / qt2d);
        if (globalMinExact > 0)
            Console.WriteLine("q_t^2 / globalMinExact           = {0} (rem {1})", qt2 / globalMinExact, qt2 % globalMinExact);
        Console.WriteLine("Minimizers found ({0}):", globalMinimizers.Count);
        foreach (var mr in globalMinimizers)
            Console.WriteLine("  g0={0:X16} g1={1:X16} g2={2:X16} g3={3:X16} spanned={4}",
                mr.g0, mr.g1, mr.g2, mr.g3, IsSpanned4(mr.g0, mr.g1, mr.g2, mr.g3));
    }
}
