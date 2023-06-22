var data = {
	"qc_metrics": {
		"n_total_pairs": 5000,
		"n_total_pairs_success": 4000,
		"n_total_pairs_failure": 1000,
		"n_total_umi": 10000,
		"n_total_cells": 375,
		"n_pairs_persample": [{
				"sample_name": "sx11_s1",
				"n_pairs_success": 350
			},
			{
				"sample_name": "sx11_s2",
				"n_pairs_success": 400
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			},
			{
				"sample_name": "sx11_s3",
				"n_pairs_success": 600
			}
		]
	}
}

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
	document.getElementById("n_total_pairs_success_perc").innerHTML = data.qc_metrics.n_total_pairs_success / data.qc_metrics.n_total_pairs * 100 + "%";
	document.getElementById("n_total_pairs_failure_perc").innerHTML = data.qc_metrics.n_total_pairs_failure / data.qc_metrics.n_total_pairs * 100 + "%";
	document.getElementById("n_total_umi").innerHTML = data.qc_metrics.n_total_umi;
	document.getElementById("n_total_cells").innerHTML = data.qc_metrics.n_total_cells;
});

