# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import functools
import logging
import time

from rinpy.constants import STAR_PRINT_COUNT


def log_time(message=None):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            msg_prefix = f"{message} " if message else "x"
            logging.info(f"{msg_prefix} started at: {time.ctime(start_time)}")
            result = func(*args, **kwargs)
            end_time = time.time()
            logging.info(f"{msg_prefix} ended at: {time.ctime(end_time)}")
            duration = end_time - start_time
            logging.info(f"{msg_prefix}Duration: {duration:.2f} seconds")
            return result

        return wrapper

    return decorator


def log_details(message):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            msg = message(*args, **kwargs) if callable(message) else message
            logging.info(f"{msg} STARTS...")
            result = func(*args, **kwargs)
            logging.info(f"{msg} ENDS...")
            return result

        return wrapper

    return decorator


def log_star_message(message: str, star_count: int = STAR_PRINT_COUNT):
    stars = "*" * star_count
    logging.info(f'{stars} {message} {stars}')


def log_with_stars(label: str):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            log_star_message(f"{label.upper()} STARTS")
            result = func(*args, **kwargs)
            log_star_message(f"{label.upper()} ENDS")
            return result

        return wrapper

    return decorator


def log_elapsed_time1(msg, start, end):
    duration = end - start
    logging.info(f"{msg} is {duration:.2f} seconds.")


def log_elapsed_time(message: str = ""):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self_obj = args[0]  # first argument to a method is always 'self'
            pdb_name = getattr(self_obj, "pdb_name", "UnknownPDB")
            start_time = time.time()
            logging.info(f'{message} STARTs at {time.ctime(start_time)}')
            result = func(*args, **kwargs)
            end_time = time.time()
            logging.info(f'{message} ENDs at {time.ctime(end_time)}')
            duration = end_time - start_time
            logging.info(f"Elapsed duration of {message} for {pdb_name.upper()} is {duration:.2f} seconds.")
            return result

        return wrapper

    return decorator
