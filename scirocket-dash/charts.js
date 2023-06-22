
document.addEventListener("DOMContentLoaded", function () { 

  window.ApexCharts &&
    new ApexCharts(document.getElementById("chart-active-users"), {
      chart: {
        type: "bar",
        fontFamily: "inherit",
        height: 40.0,
        sparkline: {
          enabled: true,
        },
        animations: {
          enabled: false,
        },
      },
      plotOptions: {
        bar: {
          columnWidth: "50%",
        },
      },
      dataLabels: {
        enabled: false,
      },
      fill: {
        opacity: 1,
      },
      series: [
        {
          name: "No. of read-pairs",
          data: n_pairs_persample_values,
        },
      ],
      tooltip: {
        theme: "dark",
      },
      grid: {
        strokeDashArray: 4,
      },
      xaxis: {
        labels: {
          padding: 0,
        },
        tooltip: {
          enabled: false,
        },
        axisBorder: {
          show: false,
        },
        type: "string",
      },
      yaxis: {
        labels: {
          padding: 4,
        },
      },
      labels: n_pairs_persample_labels,
      colors: [tabler.getColor("primary")],
      legend: {
        show: false,
      },
    }).render();
});

// Apex Charts

function generateData(count, yrange) {
  var i = 0;
  var series = [];
  while (i < count) {
    var x = (i + 1).toString();
    var y = Math.floor(Math.random() * (yrange.max - yrange.min + 1)) + yrange.min;

    series.push({
      x: x,
      y: y,
    });
    i++;
  }
  return series;
}

var options = {
  chart: {
    height: 400,
    width: 600,
    type: "heatmap",
  },

  plotOptions: {
    heatmap: {
      enableShades: true,
      colorScale: {
        ranges: [
          {
            from: -50,
            to: 0,
            color: "#d7e1ee",
          },
          {
            from: 0,
            to: 50,
            color: "#a4a2a8",
          },
          {
            from: 51,
            to: 100,
            color: "#991f17",
          },
        ],
      },
    },
  },
  colors: ["#FFFFFF"],
  dataLabels: {
    enabled: false,
  },
  series: [
    {
      name: "H",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },
    {
      name: "G",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },
    {
      name: "F",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },
    {
      name: "E",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },

    {
      name: "D",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },
    {
      name: "C",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },
    {
      name: "B",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },
    {
      name: "A",
      data: generateData(12, {
        min: -30,
        max: 55,
      }),
    },
  ],
};

document.addEventListener("DOMContentLoaded", function () {
  window.ApexCharts && new ApexCharts(document.getElementById("chart-heat"), options).render();
});
