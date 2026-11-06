import pygsheets
import pandas as pd
from pathlib import Path
import os


def usage_logging_enabled():
    """Usage logging is opt-in: GENEHACKMAN_USAGE_LOGGING=true in .env (and the developers' credentials file)."""
    return os.getenv("GENEHACKMAN_USAGE_LOGGING", "").strip().lower() in ("1", "true", "yes")


def update_google_sheet(pipeline_name, succeeded=True, error_file=None, is_test=False):
    """Opt-in usage log for the developers (see README, 'Usage logging'): user name, host, pipeline, time and
    outcome, plus the input YAML and the error log when a run fails.  Off unless usage_logging_enabled()."""
    if is_test or not usage_logging_enabled(): return

    creds_file = GENOMIC_DATA_DIR + "/googlesheets/google_sheets_creds.json"
    try:
        if not os.path.isfile(creds_file):
            return
        gc = pygsheets.authorize(service_file=creds_file)
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
            if error_file and os.path.isfile(error_file):
                with open(error_file, 'r') as file:
                    error_message = file.read()
            else:
                error_message = "(no batch log file; see Snakemake console or .snakemake/log/)"

        values = [user, hostname, pipeline_name, time_submitted_str, time_taken, result, input_data, error_message]
        worksheet.append_table(values=[values])
    except Exception as e:
        print(e)
        return

