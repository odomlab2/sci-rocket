// This file houses the code for all the charts on the dashboard.

//--------------------------------------------
// Chart - No. of succesfull pairs per sample.
//--------------------------------------------
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
        valueYField: "n_pairs_success",
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
    xAxis.data.setAll(data.qc_metrics.n_pairs_persample);
    series.data.setAll(data.qc_metrics.n_pairs_persample);

    // Add export menu.
    var exporting = am5plugins_exporting.Exporting.new(root, {
      menu: am5plugins_exporting.ExportingMenu.new(root, {}),
      dataSource: data.qc_metrics.n_pairs_persample,
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

    createChart_Heatmap(root, data_rt_plate_01, am5.color(0xe486a0));
  });
});

// p7 plate
document.addEventListener("DOMContentLoaded", function () {
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-well-p7");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data_rt_plate_01, am5.color(0x9ecfdc));
  });
});

// RT plate 01
document.addEventListener("DOMContentLoaded", function () {
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_01");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data_rt_plate_01, am5.color(0xfa1e44));
  });
});

// RT plate 02
document.addEventListener("DOMContentLoaded", function () {
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_02");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data_rt_plate_01, am5.color(0xffc825));
  });
});

// RT plate 03
document.addEventListener("DOMContentLoaded", function () {
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_03");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data_rt_plate_01, am5.color(0x5bb08f));
  });
});

// RT plate 04
document.addEventListener("DOMContentLoaded", function () {
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-rt_plate_04");
    root._logo.dispose();
    root.setThemes([am5themes_Animated.new(root)]);

    createChart_Heatmap(root, data_rt_plate_01, am5.color(0x01b4ee));
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

    createChart_Uncorrectable(root, data_top_uncorrectables_p5, am5.color(0xe95d5d));
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

    createChart_Uncorrectable(root, data_top_uncorrectables_p7, am5.color(0x426da8));
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

    createChart_Uncorrectable(root, data_top_uncorrectables_ligation, am5.color(0x297a18));
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

    createChart_Uncorrectable(root, data_top_uncorrectables_RT, am5.color(0xf27c6e));
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
      valueField: "value",
    })
  );

  series.columns.template.setAll({
    tooltipText: "{value}",
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
