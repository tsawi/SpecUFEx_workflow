# SpecUFEx Workflow

This is a workflow guide for using SpecuFEx, a waveform to clustered fingerprint algorithm introduced by Holtzman et al. (2018).

For more detailed information about the SpecUFEx algorithm, you can visit the [SpecUFEx GitHub repository](https://github.com/SpecUFex) and refer to the original publication by Holtzman et al. (2018).

## Workflow Steps

1. Create a New YAML File:
   - Create a new `.yaml` file in the `data/yaml` directory.
   - Set the required parameters in this YAML file according to your experiment.

2. Run Python Scripts:
   - Run the provided Python scripts in the `src/` directory in a specific order.
   - Start with `1_wfToSgramH5.py`, followed by `2_SpecUFEx.py`, and finally `3_Clustering.py`.
   - Provide the path to the previously created `.yaml` file as an argument.
   - For example:
     ```
     python3 1_wfToSgramH5.py ../data/yaml/demo/yaml
     ```

3. Save Figures:
   - The resulting figures are automatically saved in the `reports/figures` directory.

4. Experiment Adjustment:
   - If needed, adjust the parameters in the `.yaml` file to set up a new experiment.
   - After making changes, save the `.yaml` file.
   - Repeat the steps 2 and 3 to run and visualize the results of the new experiment.

Feel free to experiment with different settings, parameters, and data to explore the capabilities of SpecUFEx!

For more details and guidance, refer to the original resources and documentation mentioned above.

*Reference:*
- Holtzman, B. K., Paté, A., Paisley, J., Waldhauser, F., & Repetto, D. (2018). Machine learning reveals cyclic changes in seismic source spectra in Geysers geothermal field. Science Advances, 4, 5. [https://doi.org/10.1126/sciadv.aao2929](https://doi.org/10.1126/sciadv.aao2929)
