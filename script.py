import sys
import boto3
from boto3.session import Session
from datetime import datetime, timezone
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

os.environ['AWS_REGION'] = 'us-east-1'
profile_name = 'sso-admin'

session = Session(profile_name=profile_name)

dynamodb = session.resource('dynamodb',  region_name='us-east-1')
s3_client = session.client('s3', region_name='us-east-1')

DYNAMO_DB_TABLE_NAME = 'radar_profile'
table = dynamodb.Table(DYNAMO_DB_TABLE_NAME)

def extract_file_details(filename):
    """Extract radar name and timestamp from filename."""
    radar = filename[:4]
    timestamp_str = filename[4:19]  # CODE(4) + YEAR(4) + MONTH(2) + DAY(2) + "_" + HH(2) + MM(2) + SS(2) = 19 characters

    datetime_format = "%Y%m%d_%H%M%S"
    try:
        timestamp_dt = datetime.strptime(timestamp_str, datetime_format)
        timestamp_dt = timestamp_dt.replace(tzinfo=timezone.utc)
        timestamp_epoch = int(timestamp_dt.timestamp())
        return radar, timestamp_str, timestamp_epoch
    except ValueError as e:
        print(f"Error parsing timestamp for {filename}: {e}")
        return None, None, None

def process_line(line):
    """Process each line received from the Bash script."""
    logging.info(f"Processing line: {line.strip()}")
    filename = line.strip()
    if not filename or not filename.endswith('.gz'):
        logging.debug("Received non-gz file or empty line, skipping: %s", filename)
        return
    
    radar, timestamp_str, timestamp_epoch = extract_file_details(filename)
    if radar is None:
        logging.error("Skipping due to extraction failure: %s", filename)
        return
    
    try:
        response = table.put_item(
            Item={
                'radar': radar,
                'filename': filename,
                'timestamp_str': timestamp_str,
                'timestamp_epoch': timestamp_epoch,
                'created_at': int(datetime.now(timezone.utc).timestamp())
            }
        )
        logging.info("Uploaded: %s", filename)
    except Exception as e:
        logging.error("Error processing %s: %s", filename, e)

for line in sys.stdin:
    process_line(line)