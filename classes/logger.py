import logging as log
import os
# from classes.interactive import answer_yes


def create_logger(args):
    logger = log.getLogger("Logfile_logger")
    logger.setLevel(log.DEBUG)
    log_filename = "{}.log".format(args.name)
    # if os.path.isfile(log_filename):
    #     input_key = input("Logfile '{}' (default name) already exists."
    #                       "You can either overwrite this file, or specify the name by using "
    #                       "--name (-n) argument. Do you want to overwrite? (Y[es]/N[o]) > ".format(log_filename))
    #     if not answer_yes(input_key):
    #         exit()

    logfile_handler = log.FileHandler(log_filename, mode='w')
    logfile_format = log.Formatter('%(asctime)s:%(message)s')
    logfile_handler.setFormatter(logfile_format)
    logfile_handler.setLevel(log.DEBUG)
    stream_handler = log.StreamHandler()
    stream_level = log.INFO
    # if args.verbose:
    #     stream_level = log.DEBUG
    # if args.silent:
    #     stream_level = log.CRITICAL
    stream_handler.setLevel(stream_level)
    logger.addHandler(logfile_handler)
    logger.addHandler(stream_handler)
    return logger