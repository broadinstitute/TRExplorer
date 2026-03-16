"""Shared fixtures for Playwright integration tests."""

import http.server
import os
import threading

import pytest


# Directory containing the generated HTML files (index.html, locus.html, etc.)
WEBSITE_DIR = os.path.join(os.path.dirname(__file__), "..", "..")


class QuietHTTPHandler(http.server.SimpleHTTPRequestHandler):
    """HTTP handler that suppresses request logging."""

    def log_message(self, format, *args):
        pass


@pytest.fixture(scope="session")
def base_url():
    """Start a local HTTP server serving the generated website files."""
    handler = lambda *args, **kwargs: QuietHTTPHandler(*args, directory=WEBSITE_DIR, **kwargs)
    server = http.server.HTTPServer(("127.0.0.1", 0), handler)
    port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{port}"
    server.shutdown()


@pytest.fixture(scope="session")
def browser_context_args(browser_context_args):
    """Configure browser context with longer timeout and ignore HTTPS errors."""
    return {
        **browser_context_args,
        "ignore_https_errors": True,
    }


@pytest.fixture
def index_page(page, base_url):
    """Navigate to the index page and wait for it to be ready."""
    page.goto(f"{base_url}/index.html")
    page.wait_for_load_state("networkidle")
    return page


@pytest.fixture
def console_errors(page):
    """Collect JS errors emitted during the test."""
    errors = []
    page.on("console", lambda msg: errors.append(msg.text) if msg.type == "error" else None)
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    return errors
