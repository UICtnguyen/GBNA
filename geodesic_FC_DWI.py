import os
import csv
import pandas as pd 
import scipy.io as sio
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import linalg as LA

class GeodesicBrainNetworkAnalysis:
    def __init__(self, fc_dir, dwi_dir):
        """
        Initializes the class with the given directories for functional connectivity (FC) 
        and structural connectivity (DWI) data.

        Parameters:
        - fc_dir (str): Path to the directory containing functional connectivity .mat files.
        - dwi_dir (str): Path to the directory containing DWI .txt files.
        """
        self.fc_dir = fc_dir
        self.dwi_dir = dwi_dir
        self.fc_eigenvalues = {}
        self.fc_eigenvectors = {}
        self.all_R = {}
        self.DWI = {}

    def compute_fc_eigenvalues(self):
        """
        Computes the eigenvalues and eigenvectors for the functional connectivity matrices.
        """
        for subdir, dirs, files in os.walk(self.fc_dir):
            for file in files:
                if file.endswith(".mat"):
                    mat = sio.loadmat(os.path.join(subdir, file))
                    Z = mat['Z'][:166, :166]
                    np.fill_diagonal(Z, 0)
                    R = np.tanh(Z)
                    self.all_R[subdir] = R
                    w, v = np.linalg.eigh(R + np.eye(R.shape[0]))
                    w[w < 1e-13] = 1e-13
                    self.fc_eigenvalues[subdir], self.fc_eigenvectors[subdir] = w, v

    def compute_dwi_laplacian(self):
        """
        Computes the normalized Laplacian from the DWI data (density.txt files).
        """
        for root, dirs, files in os.walk(self.dwi_dir):
            for file in files:
                if file.endswith("density.txt"):
                    A = np.loadtxt(os.path.join(root, file))
                    D = np.diag(A.sum(0))
                    L = D - A  # regular Laplacian
                    D_ns = np.diag(A.sum(0)**-0.5)
                    L_n = D_ns @ L @ D_ns  # normalized Laplacian
                    self.DWI[root] = L_n

    def filter_dict_by_keys(self, dict1, dict2):
        """
        Filter keys and values in dict1 based on the keys present in dict2.
        """
        print(f"Filtering subjects... Common subjects in both datasets:")
        common_keys = set(dict1.keys()) & set(dict2.keys())
        print(f"Common subjects: {common_keys}")
        
        # Filter dict1 based on the common keys
        filtered_dict = {key: dict1[key] for key in common_keys}
        return filtered_dict

    def is_pos_def(self, x):
        """
        Check if a matrix is positive definite.
        """
        return np.all(np.linalg.eigvals(x) >= 0)

    def geo_dis(self, Q1, Q2):
        """
        Compute the geodesic distance between two matrices.

        Parameters:
        - Q1, Q2 (ndarray): Two square matrices.

        Returns:
        - geodesic_distance (float): The computed geodesic distance.
        """
        u, s, _ = LA.svd(Q1, full_matrices=True)
        Q = np.dot(u, np.dot(np.diag(s**(-1/2)), np.transpose(u)))
        M = Q @ Q2 @ Q
        _, s, _ = LA.svd(M, full_matrices=True)
        return np.sqrt(np.sum(np.log(s)**2))    

    def prepare_data(self):
        """
        Filters the subjects present in both the FC and DWI datasets.
        Strips the directory paths from the subject names and then filters them.
        """
        # Remove the directory paths from the FC dataset, correctly extracting subject names
        all_R = {os.path.basename(k): v for k, v in self.all_R.items()}
        self.values_FC = list(all_R.values())  # Extract arrays
        self.subjects_FC = list(all_R.keys())  # Extract subjects

        # Remove the directory paths from the DWI dataset, correctly extracting subject names
        DWI = {os.path.basename(k): v for k, v in self.DWI.items()}
        self.values_DWI = list(DWI.values())  # Extract arrays
        self.subjects_DWI = list(DWI.keys())  # Extract subjects

        # Debugging: Print out the raw subject identifiers
        print(f"Raw FC subjects: {self.subjects_FC}")
        print(f"Raw DWI subjects: {self.subjects_DWI}")

        # Ensure same participants for both FC and DWI datasets
        # Filtering subjects that are common to both FC and DWI datasets
        new_DWI = self.filter_dict_by_keys(DWI, all_R)
        new_FC = self.filter_dict_by_keys(all_R, DWI)

        # Debugging: Check the filtered data
        print(f"Filtered DWI subjects: {list(new_DWI.keys())}")
        print(f"Filtered FC subjects: {list(new_FC.keys())}")

        # Update the filtered FC and DWI data
        self.values_FC = list(new_FC.values())
        self.subjects_FC = list(new_FC.keys())
        self.values_DWI = list(new_DWI.values())
        self.subjects_DWI = list(new_DWI.keys())

        # Final debugging: Check the lengths of the filtered values
        print(f"Number of DWI values: {len(self.values_DWI)}")
        print(f"Number of FC values: {len(self.values_FC)}")



    def plot_geodesic_distance_matrix(self, distance_matrix, data_type):
        """
        Plot a heatmap of the geodesic distance matrix for FC or DWI data.

        Parameters:
        - distance_matrix (ndarray): A square matrix of distances.
        - data_type (str): Type of the data ('FC' or 'DWI') for labeling.
        """
        if distance_matrix.size == 0:
            print(f"Warning: The {data_type} distance matrix is empty. Skipping visualization.")
            return  # Exit if the distance matrix is empty

        sns.heatmap(distance_matrix, cmap='viridis')
        plt.title(f'{data_type} Geodesic Distance Matrix')
        plt.show()

    def compute_geodesic_distances(self):
        """
        Compute and visualize the geodesic distance matrices for both FC and DWI data.
        """
        # Calculate geodesic distance matrix for FC
        print("Computing geodesic distances for FC data...")
        self.compute_geodesic_distance_matrix(self.values_FC, "FC")

        # Calculate geodesic distance matrix for DWI
        print("Computing geodesic distances for DWI data...")
        self.compute_geodesic_distance_matrix(self.values_DWI, "DWI")

    def compute_geodesic_distances(self, data):
        """
        Computes the geodesic distance matrix for the provided data.
        """
        # Initialize the distance matrix
        n = len(data)
        dist_matrix = np.zeros((n, n))
        
        # Compute geodesic distances (example logic)
        for i in range(n):
            for j in range(i + 1, n):
                dist_matrix[i, j] = self.geo_dis(data[i], data[j])  # Assuming geo_dis is defined elsewhere
                dist_matrix[j, i] = dist_matrix[i, j]  # Symmetrize the matrix
        
        return dist_matrix

    def export_subject_data(self, dist_FC, dist_DWI):
        """
        Export subject names and similarity matrices for FC and DWI data.
        """
        # Export subject names to a CSV file
        subjects_FC = list(self.subjects_FC)  # List of FC subject names
        subjects_DWI = list(self.subjects_DWI)  # List of DWI subject names

        # Save participant names to CSV
        participants_df = pd.DataFrame({
            'subject_FC': subjects_FC,
            'subject_DWI': subjects_DWI
        })
        participants_df.to_csv('participant_names.csv', index=False)  # Change the file name as needed
        
        # Export FC distance matrix
        with open('FC_geodesic_matrix.txt', 'w') as fc_file:
            for row in dist_FC:
                fc_file.write(' '.join([str(val) for val in row]) + '\n')
        
        # Export DWI distance matrix
        with open('DWI_geodesic_matrix.txt', 'w') as dwi_file:
            for row in dist_DWI:
                dwi_file.write(' '.join([str(val) for val in row]) + '\n')

        print("Participant names and similarity matrices exported successfully.")

    def run_analysis(self):
        """
        Run the entire analysis pipeline: compute FC eigenvalues, DWI Laplacian, 
        filter data, and compute geodesic distances.
        """
        # Compute FC eigenvalues and DWI Laplacian
        self.compute_fc_eigenvalues()
        self.compute_dwi_laplacian()

        # Filter data based on common subjects
        self.prepare_data()

        # Compute and visualize geodesic distance matrices
        if len(self.values_FC) > 0 and len(self.values_DWI) > 0:
            print(f"Computing geodesic distances for FC data...")
            dist_FC = self.compute_geodesic_distances(self.values_FC)  # Pass FC data here
            print(f"Computing geodesic distances for DWI data...")
            dist_DWI = self.compute_geodesic_distances(self.values_DWI)  # Pass DWI data here
            
            # Now export subject names and similarity matrices
            self.export_subject_data(dist_FC, dist_DWI)
        else:
            print("Error: No FC or DWI data to process.")

