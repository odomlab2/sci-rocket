var data = {
  qc_metrics: {
    n_total_pairs: 5000,
    n_total_pairs_success: 4000,
    n_total_pairs_failure: 1000,
    n_total_umi: 10000,
    n_total_cells: 375,
    n_pairs_persample: [
      {
        sample_name: "sx11_s1",
        n_pairs_success: 350,
      },
      {
        sample_name: "sx11_s2",
        n_pairs_success: 400,
      },
      {
        sample_name: "sx11_s3",
        n_pairs_success: 600,
      },
      {
        sample_name: "sx11_s4",
        n_pairs_success: 600,
      },
      {
        sample_name: "sx11_s5",
        n_pairs_success: 600,
      },
      {
        sample_name: "sx11_s6",
        n_pairs_success: 600,
      },
      {
        sample_name: "sx11_s7",
        n_pairs_success: 600,
      },
      {
        sample_name: "sx11_s8",
        n_pairs_success: 600,
      },
    ],
  },
};

var n_pairs_persample = data.qc_metrics.n_pairs_persample;
var n_pairs_persample_labels = [];
var n_pairs_persample_values = [];

for (var i = 0; i < n_pairs_persample.length; i++) {
  n_pairs_persample_labels.push(n_pairs_persample[i].sample_name);
  n_pairs_persample_values.push(n_pairs_persample[i].n_pairs_success);
}

// Update numbers
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("n_totalsamples").innerHTML = n_pairs_persample.length;
  document.getElementById("n_total_pairs").innerHTML = data.qc_metrics.n_total_pairs;
  document.getElementById("n_total_pairs_success_perc").innerHTML = (data.qc_metrics.n_total_pairs_success / data.qc_metrics.n_total_pairs) * 100 + "%";
  document.getElementById("n_total_pairs_failure_perc").innerHTML = (data.qc_metrics.n_total_pairs_failure / data.qc_metrics.n_total_pairs) * 100 + "%";
  document.getElementById("n_total_umi").innerHTML = data.qc_metrics.n_total_umi;
  document.getElementById("n_total_cells").innerHTML = data.qc_metrics.n_total_cells;
});

// Data for uncorrectable barcodes
var data_top_uncorrectables_RT = [
  {
    barcode: "ATCGGGAC",
    frequency: 920,
  },
  {
    barcode: "GGGGGGGG",
    frequency: 500,
  },
  {
    barcode: "GGGGGGGA",
    frequency: 400,
  },
  {
    barcode: "GGGGGGGC",
    frequency: 300,
  },
  {
    barcode: "GGGGGGGT",
    frequency: 200,
  },
  {
    barcode: "GGGGGGAG",
    frequency: 95,
  },
  {
    barcode: "GGGGGGCG",
    frequency: 90,
  },
  {
    barcode: "GGGGGGTG",
    frequency: 80,
  },
  {
    barcode: "GGGGGGCA",
    frequency: 70,
  },
  {
    barcode: "GGGGGGCT",
    frequency: 60,
  },
  {
    barcode: "GGGAGGGT",
    frequency: 200,
  },
  {
    barcode: "TGGGGGGT",
    frequency: 200,
  },
  {
    barcode: "GGGCGGGT",
    frequency: 200,
  },
  {
    barcode: "GCCGGGGT",
    frequency: 200,
  },
  {
    barcode: "ACGGGGGT",
    frequency: 200,
  }
];
var data_top_uncorrectables_p5 = data_top_uncorrectables_RT
var data_top_uncorrectables_p7 = data_top_uncorrectables_RT
var data_top_uncorrectables_ligation = data_top_uncorrectables_RT


var data_sankey_barcodes = [
  // Good p5 -> Good p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: 5, id: "BZ-1" },
  { from: "Good p5", to: "Good p7", value: 5, id: "BZ-2" },
  { from: "Good p7", to: "Good ligation", value: 5, id: "BZ-3" },
  { from: "Good ligation", to: "Good RT", value: 5, id: "BZ-4" },

  // Good p5 -> Good p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Good p5", value: 5, id: "B6-1" },
  { from: "Good p5", to: "Good p7", value: 5, id: "B6-2" },
  { from: "Good p7", to: "Good ligation", value: 5, id: "B6-3" },
  { from: "Good ligation", to: "Bad RT", value: 5, id: "B6-4" },

  // Good p5 -> Good p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: 5, id: "B3-1" },
  { from: "Good p5", to: "Good p7", value: 5, id: "B3-2" },
  { from: "Good p7", to: "Bad ligation", value: 5, id: "B3-3" },
  { from: "Bad ligation", to: "Good RT", value: 5, id: "B3-4" },

  // Good p5 -> Good p7 -> Bad ligation -> Bad RT
  { from: "Total Bad reads", to: "Good p5", value: 5, id: "B4-1" },
  { from: "Good p5", to: "Good p7", value: 5, id: "B4-2" },
  { from: "Good p7", to: "Bad ligation", value: 5, id: "B4-3" },
  { from: "Bad ligation", to: "Bad RT", value: 5, id: "B4-4" },

  // Good p5 -> Bad p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: 10, id: "A9-1" },
  { from: "Good p5", to: "Bad p7", value: 10, id: "A9-2" },
  { from: "Bad p7", to: "Good ligation", value: 10, id: "A9-3" },
  { from: "Good ligation", to: "Good RT", value: 10, id: "A9-4" },

  // Good p5 -> Bad p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: 5, id: "B1-1" },
  { from: "Good p5", to: "Bad p7", value: 5, id: "B1-2" },
  { from: "Bad p7", to: "Bad ligation", value: 5, id: "B1-3" },
  { from: "Bad ligation", to: "Good RT", value: 5, id: "B1-4" },

  // Good p5 -> Bad p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Good p5", value: 5, id: "B2-1" },
  { from: "Good p5", to: "Bad p7", value: 5, id: "B2-2" },
  { from: "Bad p7", to: "Good ligation", value: 5, id: "B2-3" },
  { from: "Good ligation", to: "Bad RT", value: 5, id: "B2-4" },

  // Bad p5 -> Good p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: 10, id: "A1-1" },
  { from: "Bad p5", to: "Good p7", value: 10, id: "A1-2" },
  { from: "Good p7", to: "Good ligation", value: 10, id: "A1-3" },
  { from: "Good ligation", to: "Good RT", value: 10, id: "A1-4" },

  // Bad p5 -> Bad p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: 5, id: "A2-1" },
  { from: "Bad p5", to: "Bad p7", value: 5, id: "A2-2" },
  { from: "Bad p7", to: "Good ligation", value: 5, id: "A2-3" },
  { from: "Good ligation", to: "Good RT", value: 5, id: "A2-4" },

  // Bad p5 -> Bad p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: 5, id: "A3-1" },
  { from: "Bad p5", to: "Bad p7", value: 5, id: "A3-2" },
  { from: "Bad p7", to: "Bad ligation", value: 5, id: "A3-3" },
  { from: "Bad ligation", to: "Good RT", value: 5, id: "A3-4" },

  // Bad p5 -> Bad p7 -> Bad ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: 5, id: "A4-1" },
  { from: "Bad p5", to: "Bad p7", value: 5, id: "A4-2" },
  { from: "Bad p7", to: "Bad ligation", value: 5, id: "A4-3" },
  { from: "Bad ligation", to: "Bad RT", value: 5, id: "A4-4" },

  // Bad p5 -> Bad p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: 5, id: "A5-1" },
  { from: "Bad p5", to: "Bad p7", value: 5, id: "A5-2" },
  { from: "Bad p7", to: "Good ligation", value: 5, id: "A5-3" },
  { from: "Good ligation", to: "Bad RT", value: 5, id: "A5-4" },

  // Bad p5 -> Good p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: 5, id: "A6-1" },
  { from: "Bad p5", to: "Good p7", value: 5, id: "A6-2" },
  { from: "Good p7", to: "Bad ligation", value: 5, id: "A6-3" },
  { from: "Bad ligation", to: "Good RT", value: 5, id: "A6-4" },

  // Bad p5 -> Good p7 -> Bad ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: 5, id: "A7-1" },
  { from: "Bad p5", to: "Good p7", value: 5, id: "A7-2" },
  { from: "Good p7", to: "Bad ligation", value: 5, id: "A7-3" },
  { from: "Bad ligation", to: "Bad RT", value: 5, id: "A7-4" },

  // Bad p5 -> Good p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: 5, id: "A8-1" },
  { from: "Bad p5", to: "Good p7", value: 5, id: "A8-2" },
  { from: "Good p7", to: "Good ligation", value: 5, id: "A8-3" },
  { from: "Good ligation", to: "Bad RT", value: 5, id: "A8-4" },
];

var data_rt_plate_01 = [
  {
    row: "A",
    col: "1",
    value: 10,
  },
  {
    row: "B",
    col: "1",
    value: 10,
  },
  {
    row: "C",
    col: "1",
    value: 10,
  },
  {
    row: "D",
    col: "1",
    value: 10,
  },
  {
    row: "A",
    col: "2",
    value: 1,
  },
  {
    row: "B",
    col: "2",
    value: 4,
  },
  {
    row: "C",
    col: "2",
    value: 2,
  },
  {
    row: "D",
    col: "2",
    value: 6,
  },
];
