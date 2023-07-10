// This file houses the code for all the interactive data and charts on the dashboard.

//--------------------------------------------
// Helper functions
//--------------------------------------------

function add_label(root, xAxis, yAxis, x_text, y_text) {
  // Add labels to the axes
  yAxis.children.unshift(
    am5.Label.new(root, {
      text: y_text,
      fontSize: 10,
      textAlign: "center",
      y: am5.p50,
      rotation: -90,
      fontWeight: "bold",
    })
  );

  xAxis.children.push(
    am5.Label.new(root, {
      text: x_text,
      fontSize: 10,
      textAlign: "center",
      x: am5.p50,
      fontWeight: "bold",
    })
  );
}

//--------------------------------------------
// Insert data from qc_data.js
//--------------------------------------------

sample_names = Object.keys(data.samples_qc);

function roundToOne(num) {
  return +(Math.round(num + "e+1") + "e-1");
}

// Count the total number of estimated cells over samples.
var total_estimated_cells = 0;
for (var i = 0; i < sample_names.length; i++) {
  total_estimated_cells += data.samples_qc[sample_names[i]].solo_estimated_cells;
}


// Update numbers
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("version").innerHTML = data.version;
  document.getElementById("sequencing_run").innerHTML = data.sequencing_name;
  document.getElementById("n_totalsamples").innerHTML = Intl.NumberFormat("en-US").format(sample_names.length);
  document.getElementById("n_total_pairs").innerHTML = Intl.NumberFormat("en-US").format(data.n_pairs);
  document.getElementById("n_total_pairs_success_perc").innerHTML = roundToOne((data.n_pairs_success / data.n_pairs) * 100) + "%";
  document.getElementById("n_total_pairs_failure_perc").innerHTML = roundToOne((data.n_pairs_failure / data.n_pairs) * 100) + "%";
  document.getElementById("n_total_corrections").innerHTML = Intl.NumberFormat("en-US").format(data.n_corrected_p5 + data.n_corrected_p7 + data.n_corrected_ligation + data.n_corrected_rt);
  document.getElementById("n_total_cells").innerHTML = Intl.NumberFormat("en-US").format(total_estimated_cells);

  // Set the n_total_pairs_success_perc_bar width
  document.getElementById("n_total_pairs_success_perc_bar").style.width = (data.n_pairs_success / data.n_pairs) * 100 + "%";
});

//--------------------------------------------
// Chart - No. of correctable barcodes.
//--------------------------------------------

data_correctable_barcodes = [
  {
    barcode: "p5",
    frequency: data.n_corrected_p5,
  },
  {
    barcode: "p7",
    frequency: data.n_corrected_p7,
  },
  {
    barcode: "ligation",
    frequency: data.n_corrected_ligation,
  },
  {
    barcode: "rt",
    frequency: data.n_corrected_rt,
  },
];

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-n_corrections").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-n_corrections");
    root._logo.dispose();

    var chart = root.container.children.push(
      am5percent.PieChart.new(root, {
        layout: root.verticalLayout,
      })
    );

    // Add and configure Series.
    var series = chart.series.push(
      am5percent.PieSeries.new(root, {
        alignLabels: true,
        calculateAggregates: false,
        valueField: "frequency",
        categoryField: "barcode",
      })
    );

    // Set stroke of the slices.
    series.slices.template.setAll({
      strokeWidth: 2,
      stroke: am5.color(0xffffff),
      strokeOpacity: 1,
      cornerRadius: 2.5,
    });

    // Set data.
    series.data.setAll(data_correctable_barcodes);

    // Change size of labels.
    series.labels.template.setAll({
      fontSize: 10,
      text: "{barcode}:\n{frequency}",
    });

    // Set up adapters for variable slice radius
    series.slices.template.adapters.add("radius", function (radius, target) {
      var dataItem = target.dataItem;
      var high = series.getPrivate("valueHigh");

      if (dataItem) {
        var value = target.dataItem.get("valueWorking", 0);
        return (radius * value) / high;
      }
      return radius;
    });

    // Add export menu.
    var exporting = am5plugins_exporting.Exporting.new(root, {
      menu: am5plugins_exporting.ExportingMenu.new(root, {}),
      dataSource: data_correctable_barcodes,
    });

    // Animation.
    series.appear(1000);
    chart.appear(1000, 100);
  });
});

//--------------------------------------------
// Chart - No. of succesfull pairs per sample.
//--------------------------------------------

var sample_n_pairs_success = [];
for (var i = 0; i < sample_names.length; i++) {
  sample_n_pairs_success.push({
    sample_name: sample_names[i],
    frequency: data.samples_qc[sample_names[i]].n_pairs_success,
  });
}

// Sort the array by frequency
sample_n_pairs_success.sort(function (a, b) {
  return b.frequency - a.frequency;
});

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-n_pairs_sample").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-n_pairs_sample");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    // Add XY chart.
    var chart = root.container.children.push(
      am5xy.XYChart.new(root, {
        panX: true,
        panY: true,
        wheelX: "panX",
        wheelY: "zoomX",
        pinchZoomX: true,
      })
    );

    // Add cursor.
    var cursor = chart.set("cursor", am5xy.XYCursor.new(root, {}));
    cursor.lineY.set("visible", false);

    // Create axes.
    var xRenderer = am5xy.AxisRendererX.new(root, { minGridDistance: 30 });
    xRenderer.labels.template.setAll({
      rotation: -45,
      fontSize: 8,
      centerY: am5.p50,
      centerX: am5.p100,
      paddingRight: 15,
    });

    var xAxis = chart.xAxes.push(
      am5xy.CategoryAxis.new(root, {
        maxDeviation: 0.3,
        categoryField: "sample_name",
        renderer: xRenderer,
        tooltip: am5.Tooltip.new(root, {}),
      })
    );

    var yAxis = chart.yAxes.push(
      am5xy.ValueAxis.new(root, {
        maxDeviation: 0.3,
        renderer: am5xy.AxisRendererY.new(root, {
          strokeOpacity: 0.1,
        }),
      })
    );

    // Create series.
    var series = chart.series.push(
      am5xy.ColumnSeries.new(root, {
        name: "Series 1",
        xAxis: xAxis,
        yAxis: yAxis,
        valueYField: "frequency",
        sequencedInterpolation: true,
        categoryXField: "sample_name",
        tooltip: am5.Tooltip.new(root, {
          labelText: "{valueY}",
        }),
      })
    );

    // Colors.
    series.columns.template.setAll({ cornerRadiusTL: 5, cornerRadiusTR: 5, strokeOpacity: 0.5, stroke: "black", strokeWidth: 0.8 });
    series.columns.template.adapters.add("fill", function (fill, target) {
      return chart.get("colors").getIndex(series.columns.indexOf(target));
    });

    // Set data.
    xAxis.data.setAll(sample_n_pairs_success);
    series.data.setAll(sample_n_pairs_success);

    // Add labels.
    add_label(root, xAxis, yAxis, "Sample name", "Frequency");

    // Add export menu.
    var exporting = am5plugins_exporting.Exporting.new(root, {
      menu: am5plugins_exporting.ExportingMenu.new(root, {}),
      dataSource: sample_n_pairs_success,
    });

    // Animation.
    series.appear(1000);
    chart.appear(1000, 100);
  });
});

//--------------------------------------------
// Table - STARSolo summary.
//--------------------------------------------

// For each sample, generate a row in the table.
// The row will contain the following information:
// 1. Sample name
// 2. Total input reads
// 3. Total unique reads mapped to GeneFull
// 4. Estimated number of cells
// 5. Mean reads per cell
// 6. Total unique reads mapped to GeneFull for all cells.
// 7. Mean UMI per cell
// 8. Mean gene per cell
// 9. Sequencing saturation
function generate_starsolo_table(data) {
  document.getElementById("sample-starsolo-table").innerHTML = '<table class="table" id="sample-starsolo-table"> <thead> <tr> <th style="text-align:center"> <button class="table-sort" data-sort="sort-sample">Sample</button> </th> <th style="text-align:center"> <button class="table-sort" data-sort="sort-totalreads"># Reads</button> </th> <th style="text-align:center"> <button class="table-sort" data-sort="sort-totalreadsuniquegenes"># uniquely-mapped reads<br>(genes)</button> </th> <th style="text-align:center"> <button class="table-sort" data-sort="sort-totalcellsfiltered"># Cells<br>(filtered)</button> </th> <th style="text-align:center"> <button class="table-sort" data-sortß="sort-meanreadscell">Mean reads<br>(per cell)</button> </th> <th style="text-align:center"> <button class="table-sort" data-sort="sort-totalreadsuniquegenescells"># uniquely-mapped reads<br>(genes + cell)</button> </th> <th style="text-align:center"> <button class="table-sort" data-sort="sort-meanumicell">Mean UMI<br>(per cell)</button> </th> <th style="text-align:center"> <button class="table-sort" data-sort="sort-meangenescell">Mean genes<br>(per cell)</button> </th> <th style="text-align:center"> <button class="table-sort" data-sort="sort-sequencingsaturation">Sequencing<br>saturation</button> </th> </tr> </thead> <tbody class="table-tbody"> <tr></tr> </tbody> </table>'
  
  var table = document.getElementById("sample-starsolo-table");
  for (var sample in data) {
    var row = table.insertRow(-1);
    row.insertCell(0).innerHTML = sample;
    row.insertCell(1).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_reads);
    row.insertCell(2).innerHTML = Intl.NumberFormat("en-US").format(data[sample].solo_uniq_gene_reads);
    row.insertCell(3).innerHTML = Intl.NumberFormat("en-US").format(data[sample].solo_estimated_cells);
    row.insertCell(4).innerHTML = Intl.NumberFormat("en-US").format(data[sample].solo_mean_reads_per_cell);
    row.insertCell(5).innerHTML = Intl.NumberFormat("en-US").format(data[sample].solo_uniq_gene_reads_in_cells);
    row.insertCell(6).innerHTML = Intl.NumberFormat("en-US").format(data[sample].solo_mean_umi_per_cell);
    row.insertCell(7).innerHTML = Intl.NumberFormat("en-US").format(data[sample].solo_mean_gene_per_cell);

    // Add a progress bar for the solo_sequencing_saturation
    var cell = row.insertCell(8);
    var progress = document.createElement("div");
    progress.className = "row align-items-center";
    var col1 = document.createElement("div");
    col1.className = "col-12 col-lg-auto";
    col1.innerHTML = Math.round(data[sample].solo_sequencing_saturation * 100) + "%";
    var col2 = document.createElement("div");
    col2.className = "col";
    var progress_bar = document.createElement("div");
    progress_bar.className = "progress";
    progress_bar.style = "width: 3rem";
    var progress_bar_inner = document.createElement("div");
    progress_bar_inner.className = "progress-bar";
    progress_bar_inner.style = "width: " + data[sample].solo_sequencing_saturation * 100 + "%";
    progress_bar_inner.setAttribute("role", "progressbar");
    progress_bar_inner.setAttribute("aria-valuenow", data[sample].solo_sequencing_saturation);
    progress_bar_inner.setAttribute("aria-valuemin", "0");
    progress_bar_inner.setAttribute("aria-valuemax", "100");
    var span = document.createElement("span");
    span.className = "visually-hidden";
    span.innerHTML = data[sample].solo_sequencing_saturation + "%";
    progress_bar_inner.appendChild(span);
    progress_bar.appendChild(progress_bar_inner);
    col2.appendChild(progress_bar);
    progress.appendChild(col1);
    progress.appendChild(col2);
    cell.appendChild(progress);

    // Set the sorting classes of td elements
    row.cells[0].className = "sort-sample";
    row.cells[1].className = "sort-totalreads";
    row.cells[2].className = "sort-totalreadsuniquegenes";
    row.cells[3].className = "sort-totalcellsfiltered";
    row.cells[4].className = "sort-meanreadscell";
    row.cells[5].className = "sort-totalreadsuniquegenescells";
    row.cells[6].className = "sort-meanumicell";
    row.cells[7].className = "sort-meangenescell";
    row.cells[8].className = "sort-sequencingsaturation";

  }
}
document.addEventListener("DOMContentLoaded", function () {

  generate_starsolo_table(data.samples_qc);

  const list = new List("sample-starsolo-table", {
    sortClass: "table-sort",
    listClass: "table-tbody",
    valueNames: ["sort-sample", "sort-totalreads", "sort-totalreadsuniquegenes", "sort-totalcellsfiltered", "sort-meanreadscell", "sort-totalreadsuniquegenescells", "sort-meanumicell", "sort-meangenescell", "sort-sequencingsaturation"],
  });

});



//--------------------------------------------
// Table - STAR summary.
//--------------------------------------------

// For each sample, generate a row in the table.
// The row will contain the following information:
// 1. Sample name
// 2. Total input reads
// 3. Total unique reads
// 4. Total multimapping reads
// 5. Total reads filtered out.
// 6. Total chimeric reads
function generate_star_table(data) {
  document.getElementById("sample-alignment-table").innerHTML = '<table class="table" id="sample-alignment-table"> <thead> <tr> <th> <button class="table-sort" data-sort="sort-sample">Sample</button> </th> <th> <button class="table-sort" data-sort="sort-totalreads">Total reads</button> </th> <th> <button class="table-sort" data-sort="sort-uniquereads">Uniquely-mapped</button> </th> <th> <button class="table-sort" data-sort="sort-readmulti">Multimapping</button> </th> <th style="color: #d63939"> <button class="table-sort" data-sort="sort-readfilter">Filtered</button> </th> <th style="color: #d63939"> <button class="table-sort" data-sort="sort-readchimeric">Chimeric reads</button> </th> </tr> </thead> <tbody class="table-tbody"> <tr></tr> </tbody> </table>'
  
  var table = document.getElementById("sample-alignment-table");
  for (var sample in data) {
    var row = table.insertRow(-1);
    row.insertCell(0).innerHTML = sample;
    row.insertCell(1).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_reads);

    // Add a progress bar for the uniq_reads
    var cell = row.insertCell(2);
    var progress = document.createElement("div");
    progress.className = "row align-items-center";
    var col1 = document.createElement("div");
    col1.className = "col-12 col-lg-auto";
    col1.innerHTML = Math.round((data[sample].uniq_reads / data[sample].total_reads) * 100) + "%";
    var col2 = document.createElement("div");
    col2.className = "col";
    var progress_bar = document.createElement("div");
    progress_bar.className = "progress";
    progress_bar.style = "width: 3rem";
    var progress_bar_inner = document.createElement("div");
    progress_bar_inner.className = "progress-bar";
    progress_bar_inner.style = "width: " + (data[sample].uniq_reads / data[sample].total_reads) * 100 + "%";
    progress_bar_inner.setAttribute("role", "progressbar");
    progress_bar_inner.setAttribute("aria-valuenow", data[sample].uniq_reads / data[sample].total_reads);
    progress_bar_inner.setAttribute("aria-valuemin", "0");
    progress_bar_inner.setAttribute("aria-valuemax", "100");
    
    var span = document.createElement("span");
    span.className = "visually-hidden";
    span.innerHTML = (data[sample].uniq_reads / data[sample].total_reads) * 100 + "%";
    progress_bar_inner.appendChild(span);
    progress_bar.appendChild(progress_bar_inner);
    col2.appendChild(progress_bar);
    progress.appendChild(col1);
    progress.appendChild(col2);
    cell.appendChild(progress);

    // Add the other cells.
    row.insertCell(3).innerHTML = Intl.NumberFormat("en-US").format(data[sample].multimapped_reads);
    row.insertCell(4).innerHTML = Intl.NumberFormat("en-US").format(data[sample].filtered_reads);
    row.insertCell(5).innerHTML = Intl.NumberFormat("en-US").format(data[sample].chimeric_reads);

    // Set the sorting classes of td elements
    row.cells[0].className = "sort-sample";
    row.cells[1].className = "sort-totalreads";
    row.cells[2].className = "sort-uniquereads";
    row.cells[3].className = "sort-readmulti";
    row.cells[4].className = "sort-readfilter";
    row.cells[5].className = "sort-readchimeric";
  }
}
document.addEventListener("DOMContentLoaded", function () {

  generate_star_table(data.samples_qc);

  const list = new List("sample-alignment-table", {
    sortClass: "table-sort",
    listClass: "table-tbody",
    valueNames: ["sort-sample", "sort-totalreads", "sort-uniquereads", "sort-readmulti", "sort-readfilter", "sort-readchimeric"],
  });

});


//--------------------------------------------
// Chart - Ligation usage.
//--------------------------------------------

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-ligation").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-ligation");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    // Add XY chart.
    var chart = root.container.children.push(
      am5xy.XYChart.new(root, {
        panX: true,
        panY: true,
        wheelX: "panX",
        wheelY: "zoomX",
        pinchZoomX: true,
      })
    );

    // Add cursor.
    var cursor = chart.set("cursor", am5xy.XYCursor.new(root, {}));
    cursor.lineY.set("visible", false);

    // Create axes.
    var xRenderer = am5xy.AxisRendererX.new(root, { minGridDistance: 5 });
    xRenderer.labels.template.setAll({
      rotation: -90,
      centerY: am5.p50,
      centerX: am5.p100,
      paddingRight: 15,
      fontSize: 8,
    });

    var xAxis = chart.xAxes.push(
      am5xy.CategoryAxis.new(root, {
        maxDeviation: 0.3,
        categoryField: "barcode",
        renderer: xRenderer,
        tooltip: am5.Tooltip.new(root, {}),
      })
    );

    var yAxis = chart.yAxes.push(
      am5xy.ValueAxis.new(root, {
        maxDeviation: 0.3,
        renderer: am5xy.AxisRendererY.new(root, {
          strokeOpacity: 0.1,
        }),
      })
    );

    // Create series.
    var series = chart.series.push(
      am5xy.ColumnSeries.new(root, {
        name: "Series 1",
        xAxis: xAxis,
        yAxis: yAxis,
        valueYField: "frequency",
        sequencedInterpolation: true,
        categoryXField: "barcode",
        tooltip: am5.Tooltip.new(root, {
          labelText: "{valueY}",
        }),
      })
    );

    // Colors.
    series.columns.template.setAll({ cornerRadiusTL: 5, cornerRadiusTR: 5, strokeOpacity: 0.5, stroke: "black", strokeWidth: 0.8 });
    series.columns.template.adapters.add("fill", function (fill, target) {
      return chart.get("colors").getIndex(series.columns.indexOf(target));
    });

    // Set data.
    xAxis.data.setAll(data.ligation_barcode_counts);
    series.data.setAll(data.ligation_barcode_counts);

    // Add labels.
    add_label(root, xAxis, yAxis, "Ligation barcode", "Frequency");

    // Add export menu.
    var exporting = am5plugins_exporting.Exporting.new(root, {
      menu: am5plugins_exporting.ExportingMenu.new(root, {}),
      dataSource: data.ligation_barcode_counts,
    });

    // Animation.
    series.appear(1000);
    chart.appear(1000, 100);
  });
});

//--------------------------------------------
// Chart - Sankey diagram of barcodes.
//--------------------------------------------

value_TTTT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, True, True, True)";
})[0].value;
value_TTTF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, True, True, False)";
})[0].value;
value_TTFT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, True, False, True)";
})[0].value;
value_TTFF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, True, False, False)";
})[0].value;
value_TFTT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, False, True, True)";
})[0].value;
value_TFTF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, False, True, False)";
})[0].value;
value_TFFT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, False, False, True)";
})[0].value;
value_TFFF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(True, False, False, False)";
})[0].value;
value_FTTT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, True, True, True)";
})[0].value;
value_FTTF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, True, True, False)";
})[0].value;
value_FTFT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, True, False, True)";
})[0].value;
value_FTFF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, True, False, False)";
})[0].value;
value_FFTT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, False, True, True)";
})[0].value;
value_FFTF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, False, True, False)";
})[0].value;
value_FFFT = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, False, False, True)";
})[0].value;
value_FFFF = data.uncorrectables_sankey.filter(function (d) {
  return d.source == "(False, False, False, False)";
})[0].value;

var data_sankey_barcodes = [
  // Good p5 -> Good p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: value_TTTT, id: "BZ-1" },
  { from: "Good p5", to: "Good p7", value: value_TTTT, id: "BZ-2" },
  { from: "Good p7", to: "Good ligation", value: value_TTTT, id: "BZ-3" },
  { from: "Good ligation", to: "Good RT", value: value_TTTT, id: "BZ-4" },
  { from: "Good RT", to: " ", value: value_TTTT, id: "BZ-5" },

  // Good p5 -> Good p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Good p5", value: value_TTTF, id: "B1-1" },
  { from: "Good p5", to: "Good p7", value: value_TTTF, id: "B1-2" },
  { from: "Good p7", to: "Good ligation", value: value_TTTF, id: "B1-3" },
  { from: "Good ligation", to: "Bad RT", value: value_TTTF, id: "B1-4" },
  { from: "Bad RT", to: " ", value: value_TTTF, id: "B1-5" },

  // Good p5 -> Good p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: value_TTFT, id: "B3-1" },
  { from: "Good p5", to: "Good p7", value: value_TTFT, id: "B3-2" },
  { from: "Good p7", to: "Bad ligation", value: value_TTFT, id: "B3-3" },
  { from: "Bad ligation", to: "Good RT", value: value_TTFT, id: "B3-4" },
  { from: "Good RT", to: " ", value: value_TTFT, id: "B3-5" },

  // Good p5 -> Good p7 -> Bad ligation -> Bad RT
  { from: "Total Bad reads", to: "Good p5", value: value_TTFF, id: "B4-1" },
  { from: "Good p5", to: "Good p7", value: value_TTFF, id: "B4-2" },
  { from: "Good p7", to: "Bad ligation", value: value_TTFF, id: "B4-3" },
  { from: "Bad ligation", to: "Bad RT", value: value_TTFF, id: "B4-4" },
  { from: "Bad RT", to: " ", value: value_TTFF, id: "B4-5" },

  // Good p5 -> Bad p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: value_TFTT, id: "B5-1" },
  { from: "Good p5", to: "Bad p7", value: value_TFTT, id: "B5-2" },
  { from: "Bad p7", to: "Good ligation", value: value_TFTT, id: "B5-3" },
  { from: "Good ligation", to: "Good RT", value: value_TFTT, id: "B5-4" },
  { from: "Good RT", to: " ", value: value_TFTT, id: "B5-5" },

  // Good p5 -> Bad p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Good p5", value: value_TFTF, id: "B6-1" },
  { from: "Good p5", to: "Bad p7", value: value_TFTF, id: "B6-2" },
  { from: "Bad p7", to: "Bad ligation", value: value_TFTF, id: "B6-3" },
  { from: "Bad ligation", to: "Good RT", value: value_TFTF, id: "B6-4" },
  { from: "Good RT", to: " ", value: value_TFTF, id: "B6-5" },

  // Good p5 -> Bad p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Good p5", value: value_TFFT, id: "B7-1" },
  { from: "Good p5", to: "Bad p7", value: value_TFFT, id: "B7-2" },
  { from: "Bad p7", to: "Good ligation", value: value_TFFT, id: "B7-3" },
  { from: "Good ligation", to: "Bad RT", value: value_TFFT, id: "B7-4" },
  { from: "Bad RT", to: " ", value: value_TFFT, id: "B7-5" },

  // Bad p5 -> Good p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: value_TFFF, id: "A1-1" },
  { from: "Bad p5", to: "Good p7", value: value_TFFF, id: "A1-2" },
  { from: "Good p7", to: "Good ligation", value: value_TFFF, id: "A1-3" },
  { from: "Good ligation", to: "Good RT", value: value_TFFF, id: "A1-4" },
  { from: "Good RT", to: " ", value: value_TFFF, id: "A1-5" },

  // Bad p5 -> Bad p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FFTT, id: "A2-1" },
  { from: "Bad p5", to: "Bad p7", value: value_FFTT, id: "A2-2" },
  { from: "Bad p7", to: "Good ligation", value: value_FFTT, id: "A2-3" },
  { from: "Good ligation", to: "Good RT", value: value_FFTT, id: "A2-4" },
  { from: "Good RT", to: " ", value: value_FFTT, id: "A2-5" },

  // Bad p5 -> Bad p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FFFT, id: "A3-1" },
  { from: "Bad p5", to: "Bad p7", value: value_FFFT, id: "A3-2" },
  { from: "Bad p7", to: "Bad ligation", value: value_FFFT, id: "A3-3" },
  { from: "Bad ligation", to: "Good RT", value: value_FFFT, id: "A3-4" },
  { from: "Good RT", to: " ", value: value_FFFT, id: "A3-5" },

  // Bad p5 -> Bad p7 -> Bad ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FFFF, id: "A4-1" },
  { from: "Bad p5", to: "Bad p7", value: value_FFFF, id: "A4-2" },
  { from: "Bad p7", to: "Bad ligation", value: value_FFFF, id: "A4-3" },
  { from: "Bad ligation", to: "Bad RT", value: value_FFFF, id: "A4-4" },
  { from: "Bad RT", to: " ", value: value_FFFF, id: "A4-5" },

  // Bad p5 -> Bad p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FFTF, id: "A5-1" },
  { from: "Bad p5", to: "Bad p7", value: value_FFTF, id: "A5-2" },
  { from: "Bad p7", to: "Good ligation", value: value_FFTF, id: "A5-3" },
  { from: "Good ligation", to: "Bad RT", value: value_FFTF, id: "A5-4" },
  { from: "Bad RT", to: " ", value: value_FFTF, id: "A5-5" },

  // Bad p5 -> Good p7 -> Bad ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FTFT, id: "A6-1" },
  { from: "Bad p5", to: "Good p7", value: value_FTFT, id: "A6-2" },
  { from: "Good p7", to: "Bad ligation", value: value_FTFT, id: "A6-3" },
  { from: "Bad ligation", to: "Good RT", value: value_FTFT, id: "A6-4" },
  { from: "Good RT", to: " ", value: value_FTFT, id: "A6-5" },

  // Bad p5 -> Good p7 -> Bad ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FTFF, id: "A7-1" },
  { from: "Bad p5", to: "Good p7", value: value_FTFF, id: "A7-2" },
  { from: "Good p7", to: "Bad ligation", value: value_FTFF, id: "A7-3" },
  { from: "Bad ligation", to: "Bad RT", value: value_FTFF, id: "A7-4" },
  { from: "Bad RT", to: " ", value: value_FTFF, id: "A7-5" },

  // Bad p5 -> Good p7 -> Good ligation -> Bad RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FTTF, id: "A8-1" },
  { from: "Bad p5", to: "Good p7", value: value_FTTF, id: "A8-2" },
  { from: "Good p7", to: "Good ligation", value: value_FTTF, id: "A8-3" },
  { from: "Good ligation", to: "Bad RT", value: value_FTTF, id: "A8-4" },
  { from: "Bad RT", to: " ", value: value_FTTF, id: "A8-5" },

  // Bad p5 -> Good p7 -> Good ligation -> Good RT
  { from: "Total Bad reads", to: "Bad p5", value: value_FTTT, id: "A9-1" },
  { from: "Bad p5", to: "Good p7", value: value_FTTT, id: "A9-2" },
  { from: "Good p7", to: "Good ligation", value: value_FTTT, id: "A9-3" },
  { from: "Good ligation", to: "Good RT", value: value_FTTT, id: "A9-4" },
  { from: "Good RT", to: " ", value: value_FTTT, id: "A9-5" },
];

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-sankey_barcodes").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-sankey_barcodes");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    // Create series.
    var series = root.container.children.push(
      am5flow.Sankey.new(root, {
        sourceIdField: "from",
        targetIdField: "to",
        valueField: "value",
        paddingRight: 75,
        nodePadding: 20,
        idField: "id",
      })
    );

    // Add gradient to links.
    series.links.template.setAll({
      fillStyle: "gradient",
      stroke: am5.color(0x000000),
    });

    // Theme of nodes.
    series.nodes.rectangles.template.setAll({
      fillOpacity: 0.5,
      stroke: am5.color(0x000000),
      strokeWidth: 1,
      cornerRadiusTL: 4,
      cornerRadiusTR: 4,
      cornerRadiusBL: 4,
      cornerRadiusBR: 4,
    });

    // Add the values to the nodes.
    series.nodes.labels.template.setAll({
      text: "[bold]{name}[/]\n({sumOutgoing})",
    });

    // Highlight all links with the same id beginning.
    series.links.template.events.on("pointerover", function (event) {
      var dataItem = event.target.dataItem;
      var id = dataItem.get("id").split("-")[0];

      am5.array.each(series.dataItems, function (dataItem) {
        if (dataItem.get("id").indexOf(id) != -1) {
          dataItem.get("link").hover();
        }
      });
    });

    series.links.template.events.on("pointerout", function (event) {
      am5.array.each(series.dataItems, function (dataItem) {
        dataItem.get("link").unhover();
      });
    });

    // Set data.
    series.data.setAll(data_sankey_barcodes);

    // Add export menu.
    var exporting = am5plugins_exporting.Exporting.new(root, {
      menu: am5plugins_exporting.ExportingMenu.new(root, {}),
      dataSource: data_sankey_barcodes,
    });

    // Animation.
    series.appear(1000, 100);
  });
});

//--------------------------------------------
// Chart - Heatmap of 96-well plate.
//--------------------------------------------

// p5 plate
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-well-p5").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-well-p5");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data.p5_index_counts, am5.color(0xe486a0));
  });
});

// p7 plate
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-well-p7").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-well-p7");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data.p7_index_counts, am5.color(0x9ecfdc));
  });
});

// RT plate 01
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_01").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_01");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);
    // Check if data is available.
    createChart_Heatmap(root, data.rt_barcode_counts["P01"], am5.color(0xfa1e44));
  });
});

// RT plate 02
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_02").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_02");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data.rt_barcode_counts["P02"], am5.color(0xffc825));
  });
});

// RT plate 03
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_03").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_03");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data.rt_barcode_counts["P03"], am5.color(0x5bb08f));
  });
});

// RT plate 04
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_04").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_04");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data.rt_barcode_counts["P04"], am5.color(0x01b4ee));
  });
});

//--------------------------------------------
// Chart - Top 10 uncorrectable barcodes (p5).
//--------------------------------------------
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-top_uncorrectables_p5").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-top_uncorrectables_p5");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Uncorrectable(root, data.top_uncorrectables["p5"], am5.color(0xe95d5d));
  });
});

//--------------------------------------------
// Chart - Top 10 uncorrectable barcodes (p7).
//--------------------------------------------
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-top_uncorrectables_p7").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-top_uncorrectables_p7");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Uncorrectable(root, data.top_uncorrectables["p7"], am5.color(0x426da8));
  });
});

//--------------------------------------------
// Chart - Top 10 uncorrectable barcodes (ligation).
//--------------------------------------------
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-top_uncorrectables_ligation").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-top_uncorrectables_ligation");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Uncorrectable(root, data.top_uncorrectables["ligation"], am5.color(0x297a18));
  });
});

//--------------------------------------------
// Chart - Top 10 uncorrectable barcodes (RT).
//--------------------------------------------
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-top_uncorrectables_rt").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-top_uncorrectables_rt");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Uncorrectable(root, data.top_uncorrectables["rt"], am5.color(0xf27c6e));
  });
});

//--------------------------------------------
// Function - Create Uncorrectable barcode chart.
//--------------------------------------------
function createChart_Uncorrectable(root, data, color) {
  // Add XY chart.
  var chart = root.container.children.push(
    am5xy.XYChart.new(root, {
      panX: false,
      panY: false,
      wheelX: "none",
      wheelY: "none",
    })
  );

  // Add cursor.
  var cursor = chart.set("cursor", am5xy.XYCursor.new(root, {}));
  cursor.lineX.set("visible", false);

  var yRenderer = am5xy.AxisRendererY.new(root, {
    minGridDistance: 1,
  });

  // Change font of labels.
  yRenderer.labels.template.setAll({
    fontSize: 12,
    fontFamily: "Courier New",
    fill: am5.color(0x000000),
  });

  yRenderer.grid.template.set("barcode", 1);

  var yAxis = chart.yAxes.push(
    am5xy.CategoryAxis.new(root, {
      categoryField: "barcode",
      renderer: yRenderer,
      tooltip: am5.Tooltip.new(root, { themeTags: ["axis"] }),
    })
  );

  var xRenderer = am5xy.AxisRendererX.new(root, {
    minGridDistance: 10,
  });
  xRenderer.labels.template.setAll({
    fontSize: 8,
    rotation: -90,
    centerY: am5.p50,
    centerX: am5.p100,
    paddingRight: 15,
  });

  var xAxis = chart.xAxes.push(
    am5xy.ValueAxis.new(root, {
      renderer: xRenderer,
    })
  );

  // Create series.
  var series = chart.series.push(
    am5xy.ColumnSeries.new(root, {
      sequencedInterpolation: true,
      xAxis: xAxis,
      yAxis: yAxis,
      valueXField: "frequency",
      categoryYField: "barcode",
      tooltip: am5.Tooltip.new(root, {
        pointerOrientation: "left",
        labelText: "{valueX}",
      }),
    })
  );

  series.columns.template.setAll({
    cornerRadiusTR: 5,
    cornerRadiusBR: 5,
    strokeOpacity: 0,
    stroke: "black",
    fill: color,
  });

  // Set data.
  yAxis.data.setAll(data.reverse());
  series.data.setAll(data);

  // Add labels
  add_label(root, xAxis, yAxis, "Frequency", "Barcode");

  // Add export menu.
  var exporting = am5plugins_exporting.Exporting.new(root, {
    menu: am5plugins_exporting.ExportingMenu.new(root, {}),
    dataSource: data,
  });

  // Animation.
  series.appear(1000);
  chart.appear(1000, 100);
}

//--------------------------------------------
// Function - Create heatmap of 96-well plate.
//--------------------------------------------
function createChart_Heatmap(root, data, max_color) {
  // check if there's data
  if (data == undefined || data.length == 0) {
    var modal = am5.Modal.new(root, {
      content: "Plate not used.",
    });
    modal.open();
    return;
  }

  // Create XY chart
  var chart = root.container.children.push(
    am5xy.XYChart.new(root, {
      panX: false,
      panY: false,
      wheelX: "none",
      wheelY: "none",
      layout: root.verticalLayout,
    })
  );

  // Create axes and their renderers.
  var yRenderer = am5xy.AxisRendererY.new(root, {
    visible: true,
    minGridDistance: 5,
    inversed: true,
  });

  yRenderer.grid.template.set("visible", false);

  var yAxis = chart.yAxes.push(
    am5xy.CategoryAxis.new(root, {
      maxDeviation: 0,
      renderer: yRenderer,
      categoryField: "row",
    })
  );

  var xRenderer = am5xy.AxisRendererX.new(root, {
    visible: false,
    minGridDistance: 5,
    opposite: true,
  });

  xRenderer.grid.template.set("visible", false);

  var xAxis = chart.xAxes.push(
    am5xy.CategoryAxis.new(root, {
      renderer: xRenderer,
      categoryField: "col",
    })
  );

  // Create series.
  var series = chart.series.push(
    am5xy.ColumnSeries.new(root, {
      calculateAggregates: true,
      stroke: am5.color(0xffffff),
      clustered: false,
      xAxis: xAxis,
      yAxis: yAxis,
      categoryXField: "col",
      categoryYField: "row",
      valueField: "frequency",
    })
  );

  series.columns.template.setAll({
    tooltipText: "{frequency}",
    strokeOpacity: 1,
    strokeWidth: 1,
    stroke: am5.color(0xffffff),
    width: am5.percent(100),
    height: am5.percent(100),
    cornerRadiusTL: 5,
    cornerRadiusTR: 5,
    cornerRadiusBL: 5,
    cornerRadiusBR: 5,
  });

  // Show value on hover.
  series.columns.template.events.on("pointerover", function (event) {
    var di = event.target.dataItem;
    if (di) {
      heatLegend.showValue(di.get("value", 0));
    }
  });

  // Add color legend.
  series.events.on("datavalidated", function () {
    heatLegend.set("startValue", 0);
    heatLegend.set("endValue", series.getPrivate("valueHigh"));
  });

  // Set colors.
  min_color = am5.color(0xffffff);
  max_color = max_color;

  series.set("heatRules", [
    {
      target: series.columns.template,
      min: min_color,
      max: max_color,
      dataField: "value",
      key: "fill",
    },
  ]);

  // Add heat legend.
  var heatLegend = chart.bottomAxesContainer.children.push(
    am5.HeatLegend.new(root, {
      orientation: "horizontal",
      startColor: min_color,
      endColor: max_color,
      startText: "No successful reads",
      endText: "Max. no. of successful reads",
    })
  );

  // Set data
  series.data.setAll(data);

  // Specify rows and columns.
  yAxis.data.setAll([{ row: "A" }, { row: "B" }, { row: "C" }, { row: "D" }, { row: "E" }, { row: "F" }, { row: "G" }, { row: "H" }]);
  xAxis.data.setAll([{ col: "1" }, { col: "2" }, { col: "3" }, { col: "4" }, { col: "5" }, { col: "6" }, { col: "7" }, { col: "8" }, { col: "9" }, { col: "10" }, { col: "11" }, { col: "12" }]);

  // Add export menu
  var exporting = am5plugins_exporting.Exporting.new(root, {
    menu: am5plugins_exporting.ExportingMenu.new(root, {}),
    dataSource: data,
  });

  // Animation.
  chart.appear(1000, 100);
  series.appear(1000, 100);
}
