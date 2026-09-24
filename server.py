import json
import os
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse

from solver.api import calc

ROOT = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(ROOT, "static")

FILES = {
    "/": "index.html",
    "/index.html": "index.html",
    "/style.css": "style.css",
    "/app.js": "app.js",
}
TYPES = {
    "html": "text/html; charset=utf-8",
    "css": "text/css; charset=utf-8",
    "js": "application/javascript; charset=utf-8",
}


class Handler(BaseHTTPRequestHandler):
    def _send(self, code, ctype, body):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def do_GET(self):
        u = urlparse(self.path)
        name = FILES.get(u.path)
        if name is None:
            self._send(404, "text/plain; charset=utf-8", b"Not found")
            return
        try:
            with open(os.path.join(STATIC_DIR, name), "rb") as f:
                data = f.read()
        except OSError:
            self._send(404, "text/plain; charset=utf-8", b"file not found")
            return
        ext = name.rsplit(".", 1)[-1]
        self._send(200, TYPES.get(ext, "application/octet-stream"), data)

    def do_POST(self):
        u = urlparse(self.path)
        if u.path != "/api/calc":
            self._send(404, "text/plain; charset=utf-8", b"Not found")
            return
        try:
            length = int(self.headers.get("Content-Length", 0))
            raw = self.rfile.read(length)
            req = json.loads(raw.decode("utf-8"))
        except Exception as e:
            body = json.dumps({"error": "Не удалось разобрать запрос: " + str(e)}, ensure_ascii=False).encode("utf-8")
            self._send(400, "application/json; charset=utf-8", body)
            return
        try:
            resp = calc(req)
        except Exception as e:
            resp = {"error": "Внутренняя ошибка: " + str(e)}
        body = json.dumps(resp, ensure_ascii=False).encode("utf-8")
        self._send(200, "application/json; charset=utf-8", body)

    def log_message(self, fmt, *args):
        pass


def main():
    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", "8000"))
    server = ThreadingHTTPServer((host, port), Handler)
    print(f"Матричный калькулятор запущен: http://{host}:{port}")
    print("Для остановки нажмите Ctrl+C.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()