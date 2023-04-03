# Please use relative imports!
# https://www.python.org/dev/peps/pep-0328/#rationale-for-relative-imports
import sys
from qtpy.QtWidgets import QApplication
from .main_window import Window
from cad_utils.logging import enable_logging


def main():
    # Setup default C-AD logging protocol
    enable_logging()
    # TODO: Setup entrypoint code
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    app.exec_()


if __name__ == "__main__":
    main()
