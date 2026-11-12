"""A loopback stub that answers as daimon 2.5.6 does for probe 7 (port plan §4 step 9.7).

    python3 -I -B scripts/probes/p7_daimon_stub.py PORT LOGFILE

POST /v1/agents -> 201 {"id":7,"name":"ushma"}; GET /slow sleeps 3 s (the timeout case); anything else -> 404.
Every request's raw bytes (request line, headers, body) are appended to LOGFILE, so the probe can show exactly what
sandhi sent. It binds 127.0.0.1 only and exits after 20 requests or 60 s.
"""
import socket
import sys
import time

port, log = int(sys.argv[1]), sys.argv[2]
srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
srv.bind(("127.0.0.1", port))
srv.listen(8)
srv.settimeout(60)
for _ in range(20):
    try:
        c, _addr = srv.accept()
    except socket.timeout:
        break
    c.settimeout(5)
    data = b""
    try:
        while b"\r\n\r\n" not in data:
            chunk = c.recv(4096)
            if not chunk:
                break
            data += chunk
        head, _, body = data.partition(b"\r\n\r\n")
        length = 0
        for line in head.split(b"\r\n")[1:]:
            k, _, v = line.partition(b":")
            if k.strip().lower() == b"content-length":
                length = int(v.strip())
        while len(body) < length:
            chunk = c.recv(4096)
            if not chunk:
                break
            body += chunk
        with open(log, "ab") as f:
            f.write(head + b"\r\n\r\n" + body + b"\n----\n")
        first = head.split(b"\r\n")[0]
        if first.startswith(b"POST /v1/agents "):
            payload = b'{"id":7,"name":"ushma"}'
            c.sendall(b"HTTP/1.1 201 Created\r\nContent-Type: application/json\r\nContent-Length: "
                      + str(len(payload)).encode() + b"\r\nConnection: close\r\n\r\n" + payload)
        elif first.startswith(b"GET /slow "):
            time.sleep(3)
            c.sendall(b"HTTP/1.1 200 OK\r\nContent-Length: 0\r\nConnection: close\r\n\r\n")
        else:
            c.sendall(b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n")
    except OSError:
        pass
    finally:
        c.close()
srv.close()
