import pygsheets
import pandas as pd
from pathlib import Path
import os


def update_google_sheet(pipeline_name, succeeded=True, error_file=None, is_test=False):
    if is_test: return

    try:
        gc = pygsheets.authorize(service_file=GENOMIC_DATA_DIR + "/googlesheets/google_sheets_creds.json")
        spreadsheet = gc.open('GeneHackman')
        worksheet = spreadsheet[0]

        user = os.getenv('USER')
        hostname = os.getenv('HOSTNAME')
        result = "Success" if succeeded else "Failed"
        time_taken = str(datetime.now() - start_time)
        time_submitted_str = start_time.strftime('%Y-%m-%d %H:%m')

        error_message = ""
        input_data = ""
        if not succeeded:
            with open(input_file,'r') as file:
                input_data = file.read()
            with open(error_file.rstrip(), 'r') as file:
                error_message = file.read()

        values = [user, hostname, pipeline_name, time_submitted_str, time_taken, result, input_data, error_message]
        worksheet.append_table(values=[values])
    except Exception as e:
        print(e)
        return

