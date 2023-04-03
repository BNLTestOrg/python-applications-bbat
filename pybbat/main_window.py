# Please use relative imports!
# https://www.python.org/dev/peps/pep-0328/#rationale-for-relative-imports
from qtpy.QtWidgets import QLabel
from cad_ui.general import CADMainWindow


class Window(CADMainWindow):
    def __init__(self, *args, **kwargs):
        # TODO: Setup main window
        super().__init__(*args, **kwargs)
        self.setWindowTitle("pybbat")
        self.setCentralWidget(QLabel("Hello World!"))
