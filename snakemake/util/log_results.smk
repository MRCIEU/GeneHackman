import pygsheets
import pandas as pd
from pathlib import Path
import os


def update_google_sheet(pipeline_name, succeeded=True, error_file=None, is_test=False):
    if is_test: return

    gc = pygsheets.authorize(service_file=GENOMIC_DATA_DIR + "/googlesheets/google_sheets_creds.json")
    spreadsheet = gc.open('gwaspipeline')
    worksheet = spreadsheet[0]

    user = os.getenv('USER')
    hostname = os.getenv('HOSTNAME')
    result = "Success" if succeeded else "Failed"
    time_taken = str(date_submitted - datetime.now())
    time_submitted_str = time_submitted.strftime('%Y-%m-%d %H:%m')

    with open(input_file, 'r') as file:
        input_data = file.read()

    error_message = ""
    if not succeeded:
        with open(error_file, 'r') as file:
            error_message = file.read()

    values = [user, hostname, pipeline_name, input_data, time_submitted_str, time_taken, result, error_message]
    worksheet.append_table(values=[values])

