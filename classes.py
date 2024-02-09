from getpass import getpass

import defaults as g


class QUser:
    def __init__(self, uname=None, pwd=None, has_db_access=g.QUBRICSACCESS):
        self.has_db_access = has_db_access
        self.uname = None if not has_db_access else uname
        self.pwd = None if not has_db_access else pwd

    def set_credentials(self):
        if self.has_db_access:
            self.uname = input("Username: ")
            self.pwd = getpass("Password: ")

    def get_credentials(self):
        return self.uname, self.pwd

    def set_db_access(self):
        self.has_db_access = True

    def get_db_access(self):
        return self.has_db_access


class spec_collection:
    def __init__(self):
        # Might want to use a dictionary, with key the spec name
        #  We need that anyway...
        self.collection = []

    def append(self, spec):
        self.collection.append(spec)
