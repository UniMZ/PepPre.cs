namespace UniMZ.PepPre;

using Microsoft.Data.Analysis;
using UniMZ.Core;

public static partial class PepPre
{
    public static DataFrame ReadMS(string path)
    {
        // TODO: handle other format
        // should be sorted by `ScanID`.
        return path.ReadUMZ();
    }

    public static MS1 ToMS1(DataFrame df, string id = "ScanID", string peak = "Peak")
    {
        var id_ = df[id].ToArray<int>();
        var peak_ = df[peak].ToMemoryArray<Peak<double>>();
        return (id_, peak_);
    }

    public static MS2 ToMS2(DataFrame df, string id = "ScanID", string peak = "Peak", string pre = "PrecursorScan",
        string mz = "ActivationCenter", string r = "IsolationWidth")
    {
        var id_ = df[id].ToArray<int>();
        var peak_ = df[peak].ToMemoryArray<Peak<double>>();
        var pre_ = df[pre].ToArray<int>();
        var mz_ = df[mz].ToArray<double>();
        var r_ = df[r].ToArray<double>();
        return (id_, peak_, pre_, mz_, r_);
    }

    public static (MS1 ms1, MS2 ms2) ToMS(DataFrame df, string type = "ScanType", string ms1 = "MS1",
        string ms2 = "MS2")
    {
        var ms1_ = df.Filter(df.Columns[type].ElementwiseEquals(ms1));
        var ms2_ = df.Filter(df.Columns[type].ElementwiseEquals(ms2));
        return (ToMS1(ms1_), ToMS2(ms2_));
    }
}
