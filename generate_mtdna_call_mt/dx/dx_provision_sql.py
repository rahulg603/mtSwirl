# Tool for provisioning a SQL database on DNANexus.

import pyspark
import dxpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dx-init', type=str, required=True, help='Name of SQL database used to interact with Hail objects.')

if __name__ == '__main__':
    args = parser.parse_args()
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    spark.sql(f"CREATE DATABASE {args.dx_init} LOCATION  'dnax://'")