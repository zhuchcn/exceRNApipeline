from datetime import datetime


def logger(msg):
    print(
        '[ ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' ]' + msg,
        flush=True
    )