// This file houses the code for all the interactive data and charts on the dashboard.

//--------------------------------------------
// Helper functions
//--------------------------------------------

// Helper function to create an element with class and innerHTML
const createElement = (type, className, innerHTML) => {
  const element = document.createElement(type);
  element.className = className;
  if (innerHTML) element.innerHTML = innerHTML;
  return element;
  };

//--------------------------------------------
// Function to add labels to the axes
//--------------------------------------------

const add_label = (root, xAxis, yAxis, x_text, y_text) => {
  // Create a new label for the y-axis
  const yAxisLabel = am5.Label.new(root, {
    text: y_text,
    fontSize: 10,
    textAlign: "center",
    y: am5.p50,
    rotation: -90,
    fontWeight: "bold",
  });

  // Add the label to the y-axis
  yAxis.children.unshift(yAxisLabel);

  // Create a new label for the x-axis
  const xAxisLabel = am5.Label.new(root, {
    text: x_text,
    fontSize: 10,
    textAlign: "center",
    x: am5.p50,
    fontWeight: "bold",
  });

  // Add the label to the x-axis
  xAxis.children.push(xAxisLabel);
};

//--------------------------------------------
// Insert data from qc_data.js
//--------------------------------------------

const sample_names = Object.keys(data.sample_succes);

const roundToOne = num => +(Math.round(num + "e+1") + "e-1");

// Count the total number of estimated cells over samples.
let total_estimated_cells = 0;
for (let i = 0; i < sample_names.length; i++) {
  total_estimated_cells += data.sample_succes[sample_names[i]].estimated_cells;
}

// Update numbers
document.addEventListener("DOMContentLoaded", () => {
  document.getElementById("version").innerHTML = data.version;
  document.getElementById("experiment_name").innerHTML = data.experiment_name;
  document.getElementById("n_totalsamples").innerHTML = Intl.NumberFormat("en-US").format(sample_names.length);
  document.getElementById("n_total_pairs").innerHTML = Intl.NumberFormat("en-US").format(data.n_pairs);
  document.getElementById("n_total_pairs_success_perc").innerHTML = `${roundToOne((data.n_pairs_success / data.n_pairs) * 100)}%`;
  document.getElementById("n_total_pairs_failure_perc").innerHTML = `${roundToOne((data.n_pairs_failure / data.n_pairs) * 100)}%`;
  document.getElementById("n_total_corrections").innerHTML = Intl.NumberFormat("en-US").format(data.n_corrected_p5 + data.n_corrected_p7 + data.n_corrected_ligation + data.n_corrected_rt + data.n_corrected_hashing);
  document.getElementById("n_total_cells").innerHTML = Intl.NumberFormat("en-US").format(total_estimated_cells);

  // Set the n_total_pairs_success_perc_bar width
  document.getElementById("n_total_pairs_success_perc_bar").style.width = `${(data.n_pairs_success / data.n_pairs) * 100}%`;
});

//--------------------------------------------
// Chart - No. of correctable barcodes.
//--------------------------------------------
const barcodes = ["p5", "p7", "ligation", "rt", "hashing"];

const data_correctable_barcodes = barcodes.map(barcode => ({
  barcode,
  frequency: data[`n_corrected_${barcode}`],
}));

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
    frequency: data.sample_succes[sample_names[i]].n_pairs_success,
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

function generate_starsolo_table(data) {
  const table = document.getElementById("sample-starsolo-table");
  table.innerHTML = `
    <table class="tablesorter" id="sample-starsolo-table">
      <thead>
        <tr>
          <th style="text-align:center">
            Sample
          </th>
          <th style="text-align:center">
            Total input reads
          </th>
          <th style="text-align:center">
            Total mapped reads
          </th>
          <th style="text-align:center">
            Total multimapped reads
          </th>
          <th style="text-align:center">
            Total correct reads<br>(CB, UMI, gene)
          </th>
          <th style="text-align:center">
            Total exonic reads
          </th>
          <th style="text-align:center">
            Total intronic reads
          </th>
          <th style="text-align:center">
            Total intergenic reads
          </th>
          <th style="text-align:center">
            Total mitochondrial reads
          </th>
          <th style="text-align:center">
          Total exonicAS reads
          </th>
          <th style="text-align:center">
            Total intronicAS reads
          </th>
          <th style="text-align:center">
          Mean reads per cell
          </th>
          <th style="text-align:center">
          Mean genes per cell
          </th>
          <th style="text-align:center">
          Mean UMIs per cell
          </th>
          <th style="text-align:center">
          Estimated cells
          </th>
          <th style="text-align:center">
          Sequencing saturation
          </th>
        </tr>
      </thead>
      <tbody class="table-tbody">
        <tr></tr>
      </tbody>
    </table>`;

  for (const sample in data) {
    const row = table.insertRow(-1);
    row.insertCell(0).innerHTML = sample;
    row.insertCell(1).innerHTML = Intl.NumberFormat("en-US").format(data[sample].n_pairs_success);
    row.insertCell(2).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_mapped_reads);
    row.insertCell(3).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_multimapped_reads);
    row.insertCell(4).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_correct_reads_genes);
    row.insertCell(5).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_exonic_reads);
    row.insertCell(6).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intronic_reads);
    row.insertCell(7).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intergenic_reads);
    row.insertCell(8).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_mitochondrial_reads);
    row.insertCell(9).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_exonicAS_reads);
    row.insertCell(10).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intronicAS_reads);
    row.insertCell(11).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_reads_per_cell);
    row.insertCell(12).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_genes_per_cell);
    row.insertCell(13).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_umis_per_cell);
    row.insertCell(14).innerHTML = Intl.NumberFormat("en-US").format(data[sample].estimated_cells);

    // Add a progress bar for the sequencing_saturation
    const cell = row.insertCell(15);
    const progress = createElement("div", "row align-items-center");
    const col1 = createElement("div", "col-12 col-lg-auto", `${Math.round(data[sample].sequencing_saturation * 100)}%`);
    const col2 = createElement("div", "col");
    const progress_bar = createElement("div", "progress");
    progress_bar.style.width = "3rem";
    const progress_bar_inner = createElement("div", "progress-bar");
    progress_bar_inner.style.width = `${data[sample].sequencing_saturation * 100}%`;
    progress_bar_inner.setAttribute("role", "progressbar");
    progress_bar_inner.setAttribute("aria-valuenow", data[sample].sequencing_saturation);
    progress_bar_inner.setAttribute("aria-valuemin", "0");
    progress_bar_inner.setAttribute("aria-valuemax", "100");
    const span = createElement("span", "visually-hidden", `${data[sample].sequencing_saturation}%`);
    progress_bar_inner.appendChild(span);
    progress_bar.appendChild(progress_bar_inner);
    col2.appendChild(progress_bar);
    progress.appendChild(col1);
    progress.appendChild(col2);
    cell.appendChild(progress);
  }
}

document.addEventListener("DOMContentLoaded", function () {
  generate_starsolo_table(data.sample_succes);

  // Sort using tablesorter
  $("#sample-starsolo-table").tablesorter({
    sortList: [[0, 0]],
  });
  
});

//--------------------------------------------
// Table - Hashing summary.
//--------------------------------------------

// For each hashing barcode, generate a row in the table.
// The row will contain the following information:
// 1. Experiment name
// 2. Hashing Sample
// 3. Total count (correct)
// 4. Total count (corrected)
// 5. Total count (correct - upstream)
function generate_hashing_table(data) {
  const table = document.getElementById("sample-hashing-table");
  table.innerHTML = `
    <table class="tablesorter" id="sample-hashing-table">
      <thead>
        <tr>
          <th style="text-align:center">
            Hashing barcode
          </th>
          <th style="text-align:center">
            Total count (correct)
          </th>
          <th style="text-align:center">
            Total count (corrected)
          </th>
          <th style="text-align:center">
            Total count (correct - upstream)
          </th>
        </tr>
      </thead>
      <tbody class="table-tbody">
        <tr></tr>
      </tbody>
    </table>`;
      
  for (var sample in data) {
    var row = table.insertRow(-1);
    row.insertCell(0).innerHTML = sample;
    row.insertCell(1).innerHTML = Intl.NumberFormat("en-US").format(data[sample].n_correct + data[sample].n_corrected + data[sample].n_correct_upstream);
    row.insertCell(2).innerHTML = Intl.NumberFormat("en-US").format(data[sample].n_correct);
    row.insertCell(3).innerHTML = Intl.NumberFormat("en-US").format(data[sample].n_corrected);
    row.insertCell(4).innerHTML = Intl.NumberFormat("en-US").format(data[sample].n_correct_upstream);
  }
}
document.addEventListener("DOMContentLoaded", function () {
  generate_hashing_table(data.hashing);

  $("#sample-hashing-table").tablesorter({
    sortList: [[0, 0]],
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
        panY: false,
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
      paddingRight: 1,
      fontSize: 4,
    });

    // Remove grid.
    xRenderer.grid.template.setAll({ strokeOpacity: 0 });
    
    var xAxis = chart.xAxes.push(
      am5xy.CategoryAxis.new(root, {
        maxDeviation: 0.2,
        categoryField: "barcode",
        renderer: xRenderer,
        tooltip: am5.Tooltip.new(root, {}),
      })
    );

    var yRenderer = am5xy.AxisRendererY.new(root, { minGridDistance: 30 });
    yRenderer.labels.template.setAll({
      fontSize: 6,
      textAlign: "center",
      fontWeight: "bold",
      stroke: "white",
      strokeWidth: 0.5,
    });

    var yAxis = chart.yAxes.push(
      am5xy.ValueAxis.new(root, {
        maxDeviation: 0.3,
        renderer: yRenderer,
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
        sort: false,
        tooltip: am5.Tooltip.new(root, {
          labelText: "{valueY}",
        }),
      })
    );

    // Colors.
    series.columns.template.setAll({ cornerRadiusTL: 5, cornerRadiusTR: 5, strokeOpacity: 0.1, stroke: "black", strokeWidth: 0.1 });

    // Set data.
    xAxis.data.setAll(data.ligation_barcode_counts);
    series.data.setAll(data.ligation_barcode_counts);

    // Add labels.
    add_label(root, xAxis, yAxis, "Ligation barcode", "Frequency");

    // Color the first 96 bars in red
    series.columns.template.adapters.add("fill", function (fill, target) {
      if (series.columns.indexOf(target) <= 96) {
        return am5.color("#d63939");
      } 
      if (series.columns.indexOf(target) > 96 & series.columns.indexOf(target) <= 192) {
        return am5.color("#1f77b4");
      }
      if (series.columns.indexOf(target) > 192 & series.columns.indexOf(target) <= 288) {
        return am5.color("#ff7f0e");
      }
      if (series.columns.indexOf(target) > 288 & series.columns.indexOf(target) <= 384) {
        return am5.color("#2ca02c");
      }
    });

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

  // Create series with zero being the base value.
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
      valueAxis: series.yAxis,
    })
  );

  // Set data
  series.data.setAll(data);

  // Add extra hidden zero to the end of the data.
  series.data.push({ row: "I", col: "12", frequency: 0 });

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
