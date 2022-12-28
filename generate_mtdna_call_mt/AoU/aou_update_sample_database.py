import pandas as pd
import argparse
import os, re
from datetime import datetime
from google.cloud import storage


def check_table(df, subj_col):
    tf_dupe = any(df.duplicated(subset=[subj_col]))
    if tf_dupe:
        raise ValueError('ERROR: dataframe cannot have duplicates.')


def read_existing_database(path, new_suff):
    """ Imports existing database (tsv) and produces backup.
    """
    df = pd.read_csv(path, sep='\t')
    this_path_spl = os.path.splitext(path)
    new_path = this_path_spl[0] + new_suff + this_path_spl[1]
    print(f'Backing up {path} -> {new_path}...')
    df.to_csv(new_path, sep='\t', index=False)
    return df


parser = argparse.ArgumentParser()
parser.add_argument('--database-stats', type=str, default="mt_pipeline_single_2_5_stats.tsv",
                    help="Name of file containing all completed sample-level statistics. If doesn't exist, will create newly.")
parser.add_argument('--database-paths', type=str, default='mt_pipeline_single_2_5_paths.tsv',
                    help="Name of file containing paths to all completed VCF and coverage files. If doesn't exist, will create newly.")
parser.add_argument('--database-failures', type=str, default='mt_pipeline_single_2_5_failures.tsv',
                    help="Name of file containing paths to all completed VCF and coverage files. If doesn't exist, will create newly.")
parser.add_argument('--new-stats', type=str, required=True,
                    help="Name of file with sample-level stats to add to the database.")
parser.add_argument('--new-paths', type=str, required=True,
                    help="Name of file with paths to completed VCF and coverage file to add to the database.")
parser.add_argument('--new-failures', type=str, required=True,
                    help="Failure file. Okay if file points to nothing.")

if __name__ == "__main__":
    args = parser.parse_args()
    bucket_name = re.sub('^gs://', '', os.getenv("WORKSPACE_BUCKET"))
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    new_suff = '_' + datetime.now(tz=None).strftime('%d_%b_%y_%H.%M.%S')
    
    df_stats = pd.read_csv(f"gs://{bucket_name}/{args.new_stats}", sep='\t')
    check_table(df_stats, 's')
    if storage.Blob(bucket=bucket, name=args.database_stats).exists(storage_client):
        df_stats_exist = read_existing_database(f"gs://{bucket_name}/{args.database_stats}", new_suff)
        df_out = pd.concat([df_stats_exist, df_stats], axis=0, ignore_index=True)
        print(f'Adding {str(df_stats.shape[0])} records to existing stats table. New total is {str(df_out.shape[0])}.')
        if (df_stats_exist.shape[0] + df_stats.shape[0]) != df_out.shape[0]:
            raise ValueError('ERROR: the number of records in the final stats database are not as expected.')
        check_table(df_out, 's')
        df_out.to_csv(f"gs://{bucket_name}/{args.database_stats}", sep='\t', index=False)
    else:
        print(f'Stats table generated. New total is {str(df_stats.shape[0])}.')
        df_out = df_stats
        df_out.to_csv(f"gs://{bucket_name}/{args.database_stats}", sep='\t', index=False)

    df_paths = pd.read_csv(f"gs://{bucket_name}/{args.new_paths}", sep='\t')
    check_table(df_paths, 'batch')
    if storage.Blob(bucket=bucket, name=args.database_paths).exists(storage_client):
        df_paths_exist = read_existing_database(f"gs://{bucket_name}/{args.database_paths}", new_suff)
        df_out_paths = pd.concat([df_paths_exist, df_paths], axis=0, ignore_index=True)
        print(f'Adding {str(df_paths.shape[0])} records to existing paths table. New total is {str(df_out_paths.shape[0])}.')
        if (df_paths_exist.shape[0] + df_paths.shape[0]) != df_out_paths.shape[0]:
            raise ValueError('ERROR: the number of records in the final paths database are not as expected.')
        check_table(df_out_paths, 'batch')
        df_out_paths.to_csv(f"gs://{bucket_name}/{args.database_paths}", sep='\t', index=False)
    else:
        print(f'Paths table generated. New total is {str(df_paths.shape[0])}.')
        df_out_paths = df_paths
        df_out_paths.to_csv(f"gs://{bucket_name}/{args.database_paths}", sep='\t', index=False)

    # Failures. Has a unique behavior to above:
    # - If the file is already a failure, the previous record is deleted.
    # - If there is a /new success/ that was a failure, the failure is deleted.
    # - If there is a new success that is a new failure, an error is thrown.
    # - At the end, there should be no shared records between failures and successes.
    print('Processing failures...')
    if storage.Blob(bucket=bucket, name=args.new_failures).exists(storage_client):
        df_new_fail = pd.read_csv(f"gs://{bucket_name}/{args.new_failures}", sep='\t')
        if df_new_fail.shape[0] == 0:
            print(f'No new samples failed.')
            df_new_fail = None
        else:
            check_table(df_new_fail, 's')
            print(f'New failure file loaded. {str(df_new_fail.shape[0])} samples newly failed.')

            # if there is a new success that is a new failure, throw an error
            if any(df_stats.s.isin(df_new_fail.s)):
                raise ValueError('ERROR: there should not be any successful samples that also listed as failures.')
    else:
        print(f'No new samples failed.')
        df_new_fail = None
    
    if storage.Blob(bucket=bucket, name=args.database_failures).exists(storage_client):
        df_fail_exist = read_existing_database(f"gs://{bucket_name}/{args.database_failures}", new_suff)
        orig_fail_size = df_fail_exist.shape[0]
        print(f'Failed sample database with {str(orig_fail_size)} samples loaded.')
        
        # if there is a new success that was a failure, delete the failure
        df_fail_exist = df_fail_exist[~df_fail_exist['s'].isin(df_stats['s'])]
        newly_success_rm_size = df_fail_exist.shape[0]
        if orig_fail_size > newly_success_rm_size:
            print(f'{str(newly_success_rm_size)} new completed samples previously failed, and have been removed from failure database.')
        elif orig_fail_size < newly_success_rm_size:
            raise ValueError('ERROR: new samples added to failure table somehow despite filtering out samples newly completed.')
        
        if df_new_fail is not None:
            # if there are new failures that are also old failures, delete the old versions
            df_fail_exist = df_fail_exist[~df_fail_exist['s'].isin(df_new_fail['s'])]
            if newly_success_rm_size > df_fail_exist.shape[0]:
                print(f'{str(df_fail_exist.shape[0])} newly failed samples previously failed, and have been removed from old failure database.')
            elif newly_success_rm_size < df_fail_exist.shape[0]:
                raise ValueError('ERROR: new samples added to failure table somehow despite filtering out samples newly failed.')
    else:
        print(f'No existing failed sample database found.')
        df_fail_exist = None

    if df_fail_exist is None and df_new_fail is None:
        print('No failure database and no new failures found.')
    else:
        if df_fail_exist is not None and df_new_fail is not None:
            df_fail = pd.concat([df_fail_exist, df_new_fail], axis=0)
            if (df_fail_exist.shape[0] + df_new_fail.shape[0]) != df_fail.shape[0]:
                raise ValueError('ERROR: the number of records in the final failure database are not as expected.')
        elif df_fail_exist is not None:
            df_fail = df_fail_exist
        elif df_new_fail is not None:
            df_fail = df_new_fail

        check_table(df_fail, 's')
        if any(df_out.s.isin(df_fail.s)):
            raise ValueError('ERROR: there should not be any successful samples that also failures.')
        df_fail.to_csv(f"gs://{bucket_name}/{args.database_failures}", sep='\t', index=False)