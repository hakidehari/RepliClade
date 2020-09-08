

class FileStream(object):

    def __init__(filename, filepath=None):
        self.filename = filename
        self.filepath = filepath
        self.cur_open_file = None

    
    def open_file(self):
        with open(self.filename) as f:
            self.cur_open_file = f

    def write_to_file(self, msg):
        if not self.cur_open_file:
            print("File not open for editing or not specified properly")
            return
        
        