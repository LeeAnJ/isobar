"""
"""
import sys


class FluidRP10Error:

    error_index_default = 0
    error_message_default = 'not defined'
    error_location_default = 'not defined'

    def __init__(self):
        self.index = None
        self.message = None
        self.location = None

    def initiate(self):
        self.index = FluidRP10Error.error_index_default
        self.message = FluidRP10Error.error_message_default
        self.location = FluidRP10Error.error_location_default

    def get_index(self):
        return self.index

    def get_message(self):
        return self.message

    def get_location(self):
        return self.location

    def set_location(self, location='not defined'):
        self.location = location

    def print_and_terminate(self):
        print("Critical error index   : ", self.index)
        print("               message : ", self.message)
        print("program terminated at  : ", self.location)
        sys.exit()

    def print_and_continue(self):
        print("Error index   : ", self.index)
        print("      message : ", self.message)
