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


// Function to generate a matrix chart using heatmap.js
function generateMatrixChart(data_chart, ctx, rgb_x, rgb_y, rgb_z) {

    // Convert the data to a format that heatmap.js can read.
    var data_index = [];

    // Determine max frequency of data_chart
    var max_frequency = 0;
    for (var i = 0; i < data_chart.length; i++) {
      if (data_chart[i].frequency > max_frequency) {
        max_frequency = data_chart[i].frequency;
      }
    }

    for (var i = 0; i < data_chart.length; i++) {
      data_index.push({
        x: data_chart[i].col,
        y: data_chart[i].row,
        v: data_chart[i].frequency,
        v_perc: data_chart[i].frequency / max_frequency
      });
    }

    new Chart(ctx, {
      type: 'matrix',
      data: {
        datasets: [{
          data: data_index,
          // Color alpha based on frequency compared to max frequency
          backgroundColor(context) {
            var v_perc = context.dataset.data[context.dataIndex].v_perc;
            return 'rgba(' + rgb_x + ', ' + rgb_y + ', ' + rgb_z + ', ' + v_perc + ')';
          },
          borderColor() {
            return 'rgb(0, 0, 0)';
          },
          // Round the corners of the squares
          borderRadius: 2,
          borderWidth: .5,
          width(context) {
            const a = context.chart.chartArea;
            if (!a) {
              return 0;
            }
            return (a.right - a.left) / 12 - 2;
          },
          height(context) {
            const a = context.chart.chartArea;
            if (!a) {
              return 0;
            }
            return (a.bottom - a.top) / 8 - 2;
          }
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: false,
          tooltip: {
            callbacks: {
              title() {
                return '';
              },
              label(context) {
                const v = context.dataset.data[context.dataIndex];
                return [Intl.NumberFormat("en-US").format(v.v) + ' (' + v.y + v.x + ')'];
              }
            }
          }
        },
        scales: {
          x: {
            type: 'category',
            labels: ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
            offset: true,
            ticks: {
              display: true
            },
            grid: {
              display: false
            },
          },
          y: {
            type: 'category',
            labels: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
            offset: true,
            reverse: false,
            ticks: {
              display: true
            },
            grid: {
              display: false
            },
          }
        }
      }
    }
    );
}

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
// Chart - Rescues overview (pie-chart)
//--------------------------------------------

// Generate the data for the doughnut chart.
const data_correctable_barcodes = {
  labels: ["p5", "p7", "ligation", "RT", "hashing"],
  datasets: [
    {
      data: [data.n_corrected_p5, data.n_corrected_p7, data.n_corrected_ligation, data.n_corrected_rt, data.n_corrected_hashing],
      backgroundColor: ["#d63939", "#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"],
    },
  ],
};

// Remove data which is 0
for (let i = 0; i < data_correctable_barcodes.datasets[0].data.length; i++) {
  if (data_correctable_barcodes.datasets[0].data[i] === 0) {
    data_correctable_barcodes.datasets[0].data.splice(i, 1);
    data_correctable_barcodes.labels.splice(i, 1);
    data_correctable_barcodes.datasets[0].backgroundColor.splice(i, 1);
  }
}

// Generate the doughnut chart.
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rescues").innerHTML = '';
  var ctx = document.getElementById("chart-rescues").getContext("2d");
  var chart = new Chart(ctx, {
    type: "doughnut",
    data: data_correctable_barcodes,
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          position: "right",
          labels: {
            boxWidth: 15,
            font: {
              size: 10,
            },
          },
        },
      },
    },
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
  var ctx = document.getElementById("chart-n_pairs_sample").getContext("2d");

  chart_n_pair = new Chart(ctx, {
    type: "bar",
    data: {
      labels: sample_n_pairs_success.map(function (d) {
        return d.sample_name;
      }),
      datasets: [
        {
          label: "No. of succesfull pairs",
          data: sample_n_pairs_success.map(function (d) {
            return d.frequency;
          }),
          backgroundColor: "#1f77b4",
          backgroundColor: [
            'rgba(255, 99, 132, 0.9)',
            'rgba(255, 159, 64, 0.9)',
            'rgba(255, 205, 86, 0.9)',
            'rgba(75, 192, 192, 0.9)',
            'rgba(54, 162, 235, 0.9)',
            'rgba(153, 102, 255, 0.9)',
            'rgba(201, 203, 207, 0.9)'
          ],
          borderColor: [
            'rgba(0, 0, 0)',
          ],
          borderWidth: 1,
        },
      ],
    },
    options: {
      indexAxis: 'x',
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          display: false,
        },
      },
      scales: {
        x: {
          grid: {
            display: false
          },
          ticks: {
            maxRotation: 90,
            minRotation: 90,
          }
        },
      }
    },
  }); 
});

//--------------------------------------------
// Heatmaps - Index counts
//--------------------------------------------

// Well - p5
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-well-p5").innerHTML = ''
  var ctx = document.getElementById("chart-well-p5").getContext("2d");
  generateMatrixChart(data.p5_index_counts, ctx, 31, 119, 180);
});

// Well - p7
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-well-p7").innerHTML = ''
  var ctx = document.getElementById("chart-well-p7").getContext("2d");
  generateMatrixChart(data.p7_index_counts, ctx, 255, 105, 180);
});
  
// Well - RT plate 1
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_01").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_01").getContext("2d");

  // Display "No data" if there if undefined.
  if (data.rt_barcode_counts.P01 == undefined) {
    document.getElementById("rt_plate_01").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";
  }else{
    generateMatrixChart(data.rt_barcode_counts.P01, ctx, 255, 140, 0);
  }
});

// Well - RT plate 2
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_02").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_02").getContext("2d");

  // Display "No data" if there if undefined.
  if (data.rt_barcode_counts.P02 == undefined) {
    document.getElementById("rt_plate_02").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";

  }else{
    generateMatrixChart(data.rt_barcode_counts.P02, ctx, 255, 140, 0);
  }
});

// Well - RT plate 3
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_03").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_03").getContext("2d");
    // Display "No data" if there if undefined.
  if (data.rt_barcode_counts.P03 == undefined) {
    document.getElementById("rt_plate_03").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";

  }else{
    generateMatrixChart(data.rt_barcode_counts.P03, ctx, 255, 140, 0);
  }
});

// Well - RT plate 4
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_04").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_04").getContext("2d");

  if (data.rt_barcode_counts.P04 == undefined) {
    document.getElementById("rt_plate_04").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";
  }else{
    generateMatrixChart(data.rt_barcode_counts.P04, ctx, 255, 140, 0);
  }
  
});

//--------------------------------------------
// Chart - Ligation usage.
//--------------------------------------------

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-ligation").innerHTML = ''
  var ctx = document.getElementById("chart-ligation").getContext("2d");

  // Order on name LIG01, LIG02, LIG03, LIG04, LIG05, LIG06, LIG07, LIG08, LIG09, LIG10, LIG11, LIG12, etc.
  var ligation_barcode_counts = [];
  for (var i = 0; i < data.ligation_barcode_counts.length; i++) {
    ligation_barcode_counts.push({
      barcode: data.ligation_barcode_counts[i].barcode,
      frequency: data.ligation_barcode_counts[i].frequency,
      // Remove 0 from the index if it is present.
      index: data.ligation_barcode_counts[i].barcode.replace("LIG", "").replace("LIG0", ""),
    });
  }

  // Sort the array by index
  ligation_barcode_counts.sort(function (a, b) {
    return a.index - b.index;
  });


  chart_ligation = new Chart(ctx, {
    type: "bar",
    data: {
      labels: ligation_barcode_counts.map(function (d) {
        return d.barcode;
      }),
      datasets: [
        {
          label: "Frequency",
          data: ligation_barcode_counts.map(function (d) {
            return d.frequency;
          }),
          // Color the first 96 barcodes in red, the next 96 in blue, the next 96 in green and the last 96 in yellowu sing context.
          backgroundColor(context) {
            if (context.dataIndex < 96) {
              return 'rgba(255, 99, 132, 0.9)';
            } else if (context.dataIndex < 192) {
              return 'rgba(54, 162, 235, 0.9)';
            } else if (context.dataIndex < 288) {
              return 'rgba(255, 205, 86, 0.9)';
            } else {
              return 'rgba(75, 192, 192, 0.9)';
            }
          },
          borderColor: [
            'rgba(0, 0, 0)',
          ],
          borderWidth: .1,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          display: false,
        },
      },
      scales: {
        x: {
          grid: {
            display: false
          },
          ticks: {
            maxRotation: 90,
            minRotation: 90,
          }
        },
      }
    },
  });
});

//--------------------------------------------
// Chart - Top 15 uncorrectable barcodes.
//--------------------------------------------


function generateChart_uncorrectables(id, data, color) {
  document.addEventListener("DOMContentLoaded", function () {
    document.getElementById(id).innerHTML = '';
    var ctx = document.getElementById(id).getContext("2d");

    // Order on frequency
    var chartData = [];
    for (var i = 0; i < data.length; i++) {
      chartData.push({
        barcode: data[i].barcode,
        frequency: data[i].frequency,
      });
    }

    // Sort the array by frequency
    chartData.sort(function (a, b) {
      return b.frequency - a.frequency;
    });

    new Chart(ctx, {
      type: "bar",
      data: {
        labels: chartData.map(function (d) {
          return d.barcode;
        }),
        datasets: [
          {
            label: "Frequency",
            data: chartData.map(function (d) {
              return d.frequency;
            }),
            backgroundColor: color,
            borderColor: [
              'rgba(0, 0, 0)',
            ],
            borderWidth: 1,
          },
        ],
      },
      options: {
        indexAxis: 'y',
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: {
            display: false,
          },
        },
        scales: {
          x: {
            grid: {
              display: false
            },
            ticks: {
              maxRotation: 90,
              minRotation: 90,
            }
          },
          y: {
            ticks: {
              font: {
                family: "Courier New",
                size: 12,
                autoSkip: false,
              },
            },
          },
          x: {
            ticks: {
              font: {
                family: "Courier New",
                size: 10,
              },
              maxRotation: 90,
              minRotation: 90,
            },
          }
        }
      },
    });
  });
}

// Usage example:
generateChart_uncorrectables("chart-top_uncorrectables_p5", data.top_uncorrectables.p5, color = "#d63939");
generateChart_uncorrectables("chart-top_uncorrectables_p7", data.top_uncorrectables.p7, color = "#1f77b4");
generateChart_uncorrectables("chart-top_uncorrectables_lig", data.top_uncorrectables.ligation, color = "#ff7f0e");
generateChart_uncorrectables("chart-top_uncorrectables_rt", data.top_uncorrectables.rt, color = "#2ca02c");

//--------------------------------------------
// Export data to JSON table.
//--------------------------------------------

// Print data to json-table (div)
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("json-data").innerHTML = JSON.stringify(data, undefined, 2);
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
            Input reads
          </th>
          <th style="text-align:center">
            Mapped reads<br><sub>(Genome; %)</sub>
          </th>
          <th style="text-align:center">
          Mapped reads<br><sub>(Genome - Unique; %)</sub>
          </th>
          <th style="text-align:center">
          Mapped reads<br><sub>(Genes; %)</sub>
          </th>
          <th style="text-align:center">
          Mapped reads<br><sub>(Genes - Unique; %)</sub>
          </th>
          <th style="text-align:center">
            # intronic reads
          </th>
          <th style="text-align:center">
            # intergenic reads
          </th>
          <th style="text-align:center">
            # mitochondrial reads
          </th>
          <th style="text-align:center">
          # exonic (AS) reads
          </th>
          <th style="text-align:center">
            # intronic (AS) reads
          </th>
          <th style="text-align:center">
          Reads per cell<br><sub>(Mean)</sub>
          </th>
          <th style="text-align:center">
          Genes per cell<br><sub>(Mean)</sub>
          </th>
          <th style="text-align:center">
          UMI per cell<br><sub>(Mean)</sub>
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
    row.insertCell(1).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_reads);
    row.insertCell(2).innerHTML = `${roundToOne(data[sample].perc_mapped_reads_genome * 100)}%`;
    row.insertCell(3).innerHTML = `${roundToOne(data[sample].perc_unique_reads_genome_unique * 100)}%`;
    row.insertCell(4).innerHTML = `${roundToOne(data[sample].perc_mapped_reads_gene * 100)}%`;
    row.insertCell(5).innerHTML = `${roundToOne(data[sample].perc_unique_reads_gene_unique * 100)}%`;
    row.insertCell(6).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intronic_reads);
    row.insertCell(7).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intergenic_reads);
    row.insertCell(8).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_mitochondrial_reads);
    row.insertCell(9).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_exonicAS_reads);
    row.insertCell(10).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intronicAS_reads);
    row.insertCell(11).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_reads_per_cell);
    row.insertCell(12).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_genes_per_cell);
    row.insertCell(13).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_umi_per_cell);
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
// 1. Hashing Sample
// 2. Hashing barcode
// 3. Total count (summed)
// 4. Total count (correct)
// 5. Total count (corrected)
// 6. Total count (correct - upstream)
function generate_hashing_table(data) {
  const table = document.getElementById("sample-hashing-table");
  table.innerHTML = `
    <table class="tablesorter" id="sample-hashing-table">
      <thead>
        <tr>
          <th style="text-align:center">
          Sample
          </th>
          <th style="text-align:center">
            Hashing barcode
          </th>
          <th style="text-align:center">
            Total count<br><sub>(summarized)</sub>
          </th>
          <th style="text-align:center">
            Total count<br><sub>(correct)</sub>
          </th>
          <th style="text-align:center">
            Total count<br><sub>(corrected)</sub>
          </th>
          <th style="text-align:center">
          Total count<br><sub>(correct:upstream)</sub>
          </th>
        </tr>
      </thead>
      <tbody class="table-tbody">
        <tr></tr>
      </tbody>
    </table>`;
      
  for (var sample in data) {
    for (var hashing_barcode in data[sample]) {
      const row = table.insertRow(-1);
      row.insertCell(0).innerHTML = sample;
      row.insertCell(1).innerHTML = hashing_barcode;
      row.insertCell(2).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_correct + data[sample][hashing_barcode].n_corrected + data[sample][hashing_barcode].n_correct_upstream);
      row.insertCell(3).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_correct);
      row.insertCell(4).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_corrected);
      row.insertCell(5).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_correct_upstream);
    }
  }
}

document.addEventListener("DOMContentLoaded", function () {
  generate_hashing_table(data.hashing);

  $("#sample-hashing-table").tablesorter({
    sortList: [[0, 0]],
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
  document.getElementById("chart-sankey").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-sankey");
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

    // Animation.
    series.appear(1000, 100);
  });
});
