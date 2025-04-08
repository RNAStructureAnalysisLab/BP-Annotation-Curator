# AUTHOR: Kristopher Church

# Scrapes through RNA 3D Motif Atlas for each cluster. For each cluster, the 
# web page has an HTML table, and it parses through the information in this
# table into a CSV file representation

# =============================================================================

import requests
import time
import os
import shutil
import csv
from bs4 import BeautifulSoup

class Cluster_Downloader:
    CRAWL_DELAY = 10 # https://rna.bgsu.edu/robots.txt suggests 10 seconds
    CLUSTER_DOWNLOAD_DIRECTORY = 'Data/Raw/R3DMA'
    HOMEPAGE_URL = 'https://rna.bgsu.edu/rna3dhub/motifs'
    
    failed_links = [] # R3DMA links that still need their contents downloaded
    clustered_pdb_ids = set() # PDB IDs that have been seen in downloaded data
    
    @staticmethod
    def download():
        # Ensure directory is blank
        if os.path.exists(Cluster_Downloader.CLUSTER_DOWNLOAD_DIRECTORY):
            shutil.rmtree(Cluster_Downloader.CLUSTER_DOWNLOAD_DIRECTORY)
        os.makedirs(Cluster_Downloader.CLUSTER_DOWNLOAD_DIRECTORY)
        
        # Dynamically find the URLs for all motif clusters
        urls = Cluster_Downloader.get_links(Cluster_Downloader.HOMEPAGE_URL)
        
        # Scrape motif cluster table data for each cluster from a given URL
        for url in urls:
            links = Cluster_Downloader.get_links(url)
            Cluster_Downloader.get_data_from_links(links)
 
        # Document currently used PDBs and clusters that couldn't be downloaded
        Cluster_Downloader.write_file(
            'pdb_ids.txt', Cluster_Downloader.clustered_pdb_ids
        )
        Cluster_Downloader.write_file(
            'failed_links.txt', Cluster_Downloader.failed_links
        )
        
    
    @staticmethod
    def get_links(url):
        attempts = 0
        while attempts < 3:
            response = requests.get(Cluster_Downloader.HOMEPAGE_URL)
            
            if response.status_code != 200:
                attempts += 1
            else:
                soup = BeautifulSoup(response.text, 'html.parser')
                return [
                    anchor['href'] for anchor in soup.find_all('a', href=True) 
                    if anchor['href'].startswith(
                        'https://rna.bgsu.edu/rna3dhub/motif/view/'
                    ) or anchor['href'].startswith(
                        'https://rna.bgsu.edu/rna3dhub/motifs/release/'
                    )
                ]
        
        Cluster_Downloader.failed_links.append(url)
        return []
    
    @staticmethod
    def get_data_from_links(cluster_download_links):
        session = requests.Session()
        for link in cluster_download_links:
            # Attempts to download the table for a cluster, up to 3 attempts
            for attempt in [1, 2, 3]:
                response = session.get(link)
                try:
                    if response.status_code == 200:
                        # Get header for the table
                        soup = BeautifulSoup(response.text, 'html.parser')
                        table = soup.find('table')
                        headers = [
                            table_header.text.strip() for table_header in 
                            table.find_all('th')
                        ]
                        
                        # Get rows in the table
                        rows = []
                        index = headers.index('PDB')
                        for table_row in table.find_all('tr')[1:]:
                            entries = [
                                table_data.text.strip() for table_data in 
                                table_row.find_all('td')
                            ]
                            rows.append(entries)
                            
                            # Document PDB IDs we come across for storing in a 
                            # file later
                            Cluster_Downloader.clustered_pdb_ids.add(
                                entries[index]
                            )
                        
                        # Export table contents as a CSV
                        file_name = f'{link.split('/')[-1]}.csv'
                        file_path = os.path.join(
                            Cluster_Downloader.CLUSTER_DOWNLOAD_DIRECTORY, 
                            file_name
                        )
                        with open(file_path, mode='w', newline='') as file:
                            writer = csv.writer(file)
                            
                            writer.writerow(headers)
                            writer.writerows(rows)
                        
                        time.sleep(Cluster_Downloader.CRAWL_DELAY)
                        break
                    
                    raise Exception('Could not access link')
                
                except Exception:
                    if attempt == 3: # If is final attempt
                        Cluster_Downloader.failed_links.append(link)
                    else:
                        time.sleep(Cluster_Downloader.CRAWL_DELAY)
        
    @staticmethod
    def retry_download():
        urls = Cluster_Downloader.failed_links
        Cluster_Downloader.failed_links = []
        
        # Check if the homepage is a failed link
        if Cluster_Downloader.HOMEPAGE_URL in Cluster_Downloader.failed_links:
            # Dynamically find the URLs for all motif clusters
            urls = Cluster_Downloader.get_links(
                Cluster_Downloader.HOMEPAGE_URL
            )
            
            for url in urls:
                links = Cluster_Downloader.get_links(url)
                Cluster_Downloader.get_data_from_links(links)
        else:
           Cluster_Downloader.get_data_from_links(urls)
        
        # Document currently used PDBs and clusters that couldn't be downloaded
        Cluster_Downloader.write_file(
            'pdb_ids.txt', Cluster_Downloader.clustered_pdb_ids
        )
        Cluster_Downloader.write_file(
            'failed_links.txt', Cluster_Downloader.failed_links
        )
        
    @staticmethod
    def write_file(file_name, data):
        file_path = os.path.join(
            Cluster_Downloader.CLUSTER_DOWNLOAD_DIRECTORY, file_name
        )
        with open(file_path, 'w') as file:
            for entry in data:
                file.write(entry + '\n')
        