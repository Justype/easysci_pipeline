import requests
import os

# https://genome.ucsc.edu/cgi-bin/hgTables?
# clade=mammal&
# org=Human&
# db=hg38&
# hgta_group=allTracks&
# hgta_track=rmsk&
# hgta_table=rmsk&
# hgta_regionType=genome&
# position=&
# hgta_outputType=gff&
# hgta_outFileName=hg38_rmsk.gtf&
# hgta_doTopSubmit=get+output

def main():
    # Download for Human (hg38)
    if download_ucsc_repeatmasker_gtf(genome="hg38", output_file="hg38_rmsk.gtf"):
        print("\n--- hg38 RepeatMasker GTF download completed. ---")
    else:
        print("\n--- hg38 RepeatMasker GTF download failed. ---")
    
    # Download for Mouse (mm39)
    if download_ucsc_repeatmasker_gtf(genome="mm39", output_file="mm39_rmsk.gtf"):
        print("\n--- mm39 RepeatMasker GTF download completed. ---")
    else:
        print("\n--- mm39 RepeatMasker GTF download failed. ---")

def download_ucsc_repeatmasker_gtf(genome="hg38", output_file="hg38_rmsk.gtf"):
    """
    Downloads the RepeatMasker table as a GTF file from UCSC Genome Browser.

    Args:
        genome (str): The UCSC genome assembly (e.g., "hg38", "mm10", "mm39").
        output_file (str): The name of the output GTF file.
    """
    ucsc_url = "https://genome.ucsc.edu/cgi-bin/hgTables"

    # Base parameters derived from inspecting the UCSC Table Browser
    # 'position' parameter is removed as it's not strictly necessary when regionType is 'genome'
    params = {
        'hgsid': '0', # Session ID, '0' often works.
        'clade': 'mammal', # Clade for the genome
        'org': 'Human', # Organism, can be adjusted based on genome
        'db': genome,
        'hgta_group': 'allTracks', # Group of tracks
        'hgta_track': 'rmsk', # RepeatMasker track
        'hgta_table': 'rmsk', 
        'hgta_regionType': 'genome', # Entire genome
        'hgta_outputType': 'gff', # Output format
        'hgta_outFileName': output_file, # Desired output filename
        'hgta_doTopSubmit': 'get+output', # The button click
    }

    # Adjust 'org' and 'clade' based on the specific genome assembly
    if genome == "hg38":
        params['org'] = 'Human'
        params['clade'] = 'mammal'
    elif genome == "mm39":
        params['org'] = 'Mouse'
        params['clade'] = 'mammal'
    elif genome == "mm10": # Added for completeness as it's common too
        params['org'] = 'Mouse'
        params['clade'] = 'mammal'
    # Add more 'elif' blocks for other specific genomes you might need in the future
    else:
        print(f"Warning: Genome '{genome}' is not recognized. Defaulting to Human (hg38) parameters.")
        exit(1)


    print(f"Attempting to download RepeatMasker GTF for {genome}...")
    print(f"Target URL: {ucsc_url}")
    print(f"Parameters: {params}") # Print updated parameters

    try:
        response = requests.post(ucsc_url, data=params, stream=True, timeout=120) # Added timeout
        response.raise_for_status() # Raise an exception for HTTP errors (4xx or 5xx)

        # Check for HTML response indicating an error or redirection
        if response.headers.get('Content-Type') and 'text/html' in response.headers['Content-Type']:
            print("Warning: Response content type is HTML. Download might have failed or redirected.")
            print("This usually means the parameters (especially 'db', 'org', 'clade') are incorrect or the 'hgsid' is invalid.")
            print("Partial HTML content (first 500 chars):")
            print(response.text[:500])
            return False

        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Successfully downloaded {output_file}")
        return True

    except requests.exceptions.Timeout:
        print(f"Error: The request timed out after 120 seconds for genome {genome}.")
        return False
    except requests.exceptions.RequestException as e:
        print(f"An HTTP/network error occurred: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

if __name__ == "__main__":
    main()
