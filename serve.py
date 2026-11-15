"""Serve this folder at http://localhost:8000 with caching switched off, so edits show up on a normal refresh.

Same as `python -m http.server 8000`, plus a Cache-Control header. Run from anywhere:  python serve.py
Stop with Ctrl+C. Optional: python serve.py 8080  (different port).
"""
import os
import sys
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8000


class NoCacheHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        # tell the browser never to reuse a stored copy of any file
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()


if __name__ == "__main__":
    handler = partial(NoCacheHandler, directory=HERE)       # always serve the folder this file lives in
    with ThreadingHTTPServer(("127.0.0.1", PORT), handler) as server:   # 127.0.0.1: this computer only
        print(f"Serving {HERE}\nOpen http://localhost:{PORT}  (Ctrl+C to stop)")
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nStopped.")
