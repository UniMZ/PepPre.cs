namespace UniMZ.PepPreTest;

using Microsoft.Data.Analysis;
using Microsoft.ML.Data;
using UniMZ.Core;
using UniMZ.PepPre;

[TestClass]
public class DataTest
{
    private static readonly string[] Type = ["MS1", "MS2", "MS2", "MS1", "MS2", "MS2"];
    private static readonly int[] ID = [2, 4, 6, 7, 8, 9];

    private static readonly VBuffer<Peak<double>>[] Peak =
    [
        new(2, [new Peak<double>(1.1, 1), new Peak<double>(1.2, 2)]),
        new(2, [new Peak<double>(2.1, 1), new Peak<double>(2.2, 2)]),
        new(3, [new Peak<double>(3.1, 1), new Peak<double>(3.2, 2), new Peak<double>(3.3, 3)]),
        new(2, [new Peak<double>(1.2, 1), new Peak<double>(1.3, 2)]),
        new(2, [new Peak<double>(2.1, 1), new Peak<double>(2.2, 2)]),
        new(3, [new Peak<double>(3.1, 1), new Peak<double>(3.2, 2), new Peak<double>(3.3, 3)])
    ];

    private static readonly int[] Pre = [0, 2, 2, 0, 8, 9];
    private static readonly double[] MZ = [0.0, 1.1, 1.2, 0, 1.0, 2.0];
    private static readonly double[] R = [0.0, 0.05, 0.15, 0.0, 0.05, 1.0];

    private static DataFrame BuildDataFrame()
    {
        var type = DataFrameColumn.Create("ScanType", Type);
        var id = DataFrameColumn.Create("ScanID", ID);
        var peak = new VBufferDataFrameColumn<Peak<double>>("Peak", Peak);
        var pre = DataFrameColumn.Create("PrecursorScan", Pre);
        var mz = DataFrameColumn.Create("ActivationCenter", MZ);
        var r = DataFrameColumn.Create("IsolationWidth", R);
        return new DataFrame(type, id, peak, pre, mz, r);
    }

    private static int[] GetMS1() { return Enumerable.Range(0, Type.Length).Where(i => Type[i] == "MS1").ToArray(); }
    private static int[] GetMS2() { return Enumerable.Range(0, Type.Length).Where(i => Type[i] == "MS2").ToArray(); }

    [TestMethod]
    public void TestToMS()
    {
        var df = BuildDataFrame();
        var (ms1, ms2) = PepPre.ToMS(df);
        var idx1 = GetMS1();
        var idx2 = GetMS2();
        CollectionAssert.AreEqual(ms1.id, idx1.Select(i => ID[i]).ToArray());
        CollectionAssert.AreEqual(ms2.id, idx2.Select(i => ID[i]).ToArray());
        Assert.AreEqual(idx1.Length, ms1.id.Length);
        for (var i = 0; i < idx1.Length; i++)
            CollectionAssert.AreEqual(ms1.peak[i].ToArray(), Peak[idx1[i]].GetValues().ToArray());
        Assert.AreEqual(idx2.Length, ms2.id.Length);
        for (var i = 0; i < idx2.Length; i++)
            CollectionAssert.AreEqual(ms2.peak[i].ToArray(), Peak[idx2[i]].GetValues().ToArray());
        CollectionAssert.AreEqual(ms2.pre, idx2.Select(i => Pre[i]).ToArray());
        CollectionAssert.AreEqual(ms2.mz, idx2.Select(i => MZ[i]).ToArray());
        CollectionAssert.AreEqual(ms2.r, idx2.Select(i => R[i]).ToArray());
    }
}
