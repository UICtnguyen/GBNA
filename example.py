from geodesic_FC_DWI import GeodesicBrainNetworkAnalysis

# Set paths to your data, fc = functional connectomes, dwi = structural connectomes 
fc_dir = 'your/path/here' #example: '/Users/tnguy271/Downloads/unmasck_imaging_files/rsNetwork/visit_one/'
dwi_dir = 'your/path/here' #example: '/Users/tnguy271/Downloads/unmasck_imaging_files/final_network_two/visit_one/'

# Initialize the analysis class
analysis = GeodesicBrainNetworkAnalysis(fc_dir, dwi_dir)

# Run the analysis pipeline
analysis.run_analysis()
