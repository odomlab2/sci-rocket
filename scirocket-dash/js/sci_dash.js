// This file houses the code for all the interactive data and charts on the dashboard.

//--------------------------------------------
// Insert data from qc_data.js
//--------------------------------------------

sample_names = Object.keys(data.samples_qc);

// Update numbers
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("version").innerHTML = data.version;
  document.getElementById("sequencing_run").innerHTML = data.sequencing_name;
  document.getElementById("n_totalsamples").innerHTML = Intl.NumberFormat("en-US").format(sample_names.length);
  document.getElementById("n_total_pairs").innerHTML = Intl.NumberFormat("en-US").format(data.n_pairs);
  document.getElementById("n_total_pairs_success_perc").innerHTML = (data.n_pairs_success / data.n_pairs) * 100 + "%";
  document.getElementById("n_total_pairs_failure_perc").innerHTML = (data.n_pairs_failure / data.n_pairs) * 100 + "%";
  document.getElementById("n_total_corrections").innerHTML = Intl.NumberFormat("en-US").format(data.n_corrected_p5 + data.n_corrected_p7 + data.n_corrected_ligation + data.n_corrected_rt);
  document.getElementById("n_total_cells").innerHTML = Intl.NumberFormat("en-US").format(data.n_cells);

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
      text: "{barcode}: {frequency}",
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
      rotation: -90,
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
// Table - Sample summary.
//--------------------------------------------

// For each sample, generate a row in the table.
// The row will contain the following information:
// 1. Sample name
// 2. Number of reads
// 3. Number of unique UMI
// 4. Number of cells
// 5. Number of cells with > 100 UMI
// 6. Number of cells with > 1000 UMI
// 7. Duplication rate + progress bar
function generate_sample_summary_table(data) {
  var table = document.getElementById("sample-summary-table");
  for (var sample in data) {
      var row = table.insertRow(-1);
      row.insertCell(0).innerHTML = sample;
      row.insertCell(1).innerHTML = data[sample].n_pairs_success;
      row.insertCell(2).innerHTML = data[sample].n_UMIs;
      row.insertCell(3).innerHTML = data[sample].n_cells;
      row.insertCell(4).innerHTML = data[sample].n_cells_umi_100;
      row.insertCell(5).innerHTML = data[sample].n_cells_umi_1000;

      // Add a progress bar for the duplication rate
      var cell = row.insertCell(6);
      var progress = document.createElement("div");
      progress.className = "row align-items-center";
      var col1 = document.createElement("div");
      col1.className = "col-12 col-lg-auto";
      col1.innerHTML = data[sample].duplication_rate * 100 + "%";
      var col2 = document.createElement("div");
      col2.className = "col";
      var progress_bar = document.createElement("div");
      progress_bar.className = "progress";
      progress_bar.style = "width: 5rem";
      var progress_bar_inner = document.createElement("div");
      progress_bar_inner.className = "progress-bar";
      progress_bar_inner.style = "width: " + data[sample].duplication_rate * 100 + "%";
      progress_bar_inner.setAttribute("role", "progressbar");
      progress_bar_inner.setAttribute("aria-valuenow", data[sample].duplication_rate);
      progress_bar_inner.setAttribute("aria-valuemin", "0");
      progress_bar_inner.setAttribute("aria-valuemax", "100");
      progress_bar_inner.setAttribute("aria-label", data[sample].duplication_rate + "% Complete");
      var span = document.createElement("span");
      span.className = "visually-hidden";
      span.innerHTML = data[sample].duplication_rate + "% Complete";
      progress_bar_inner.appendChild(span);
      progress_bar.appendChild(progress_bar_inner);
      col2.appendChild(progress_bar);
      progress.appendChild(col1);
      progress.appendChild(col2);
      cell.appendChild(progress);

      // Set the sorting classes of td elements
      row.cells[0].className = "sort-sample";
      row.cells[1].className = "sort-reads";
      row.cells[2].className = "sort-unique_umi";
      row.cells[3].className = "sort-totalcells";
      row.cells[4].className = "sort-totalcells_100";
      row.cells[5].className = "sort-totalcells_1000";
      row.cells[6].className = "sort-duplication";
      row.cells[6].setAttribute("data-progress", data[sample].duplication_rate);
      
  }
}
document.addEventListener("DOMContentLoaded", function () {
  generate_sample_summary_table(data.samples_qc);
});


//--------------------------------------------
// Chart - Ligation usage.
//--------------------------------------------

document.addEventListener("DOMContentLoaded", function () {
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
    var xRenderer = am5xy.AxisRendererX.new(root, { minGridDistance: 30 });
    xRenderer.labels.template.setAll({
      rotation: -90,
      centerY: am5.p50,
      centerX: am5.p100,
      paddingRight: 15,
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
document.addEventListener("DOMContentLoaded", function () {
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
    fontSize: 16,
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

  var xAxis = chart.xAxes.push(
    am5xy.ValueAxis.new(root, {
      renderer: am5xy.AxisRendererX.new(root, {
        strokeOpacity: 0.1,
      }),
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
