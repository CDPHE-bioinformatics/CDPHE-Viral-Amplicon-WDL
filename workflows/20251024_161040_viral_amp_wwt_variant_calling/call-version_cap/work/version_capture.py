import argparse
import json
import pandas as pd
import sys

def get_options(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis_date', 
                        help = 'date to add to the analysis_date column in the output CSV',
                        required = True)
    parser.add_argument('--docker_name', 
                        help = 'the Dockerfile ENV variable "NAME" for this container',
                        required = True)
    parser.add_argument('--docker_host', 
                        help = 'the Dockerfile ENV variable "HOST" for this container',
                        required = True)
    parser.add_argument('--docker_version', 
                        help = 'the Dockerfile ENV variable "VERSION" for this container',
                        required = True)
    parser.add_argument('--out_fn', 
                        help = 'the name for the output file',
                        required = True)
    parser.add_argument('--project_name', 
                        help = 'project name of batch being analyzed (e.g. cov_2022_grid)',
                        required = True)
    parser.add_argument('--versions_json', 
                        help = 'JSON file with an array of version_info objects (keys should be "software", "docker", and "version")',
                        required = True)
    parser.add_argument('--sample_name', 
                        help = 'provide only when calling from a sample-level workflow', 
                        required = False,
                        default = None)
    options = parser.parse_args(args)
    return options

def create_version_df(options):
    with open(options.versions_json) as f:
        records = json.load(f)['versions']

    docker_full_name = options.docker_host + '/' + options.docker_name

    version_capture_docker_record = {
        'software': options.docker_name,
        'docker': docker_full_name,
        'version': options.docker_version
    }
    records.insert(0, version_capture_docker_record)

    df = pd.DataFrame(records)
    df = df.rename(columns={'docker': 'associated_docker_container'})  # backwards-compat naming
    
    if options.sample_name is not None:
        insert_dict = {'sample_name': options.sample_name}
    else:
        insert_dict = {}
    insert_dict.update({'project_name': options.project_name, 
                        'analysis_date': options.analysis_date})
    for i, (k, v) in enumerate(insert_dict.items()):
        df.insert(i, k, v)
    
    return df

if __name__ == '__main__':
    options = get_options()
    df = create_version_df(options)
    df.to_csv(options.out_fn, index=False)
