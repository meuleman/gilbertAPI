import pyarrow as pa
import pyarrow.parquet as pq
from datetime import datetime, timedelta
from pathlib import Path
import threading
from typing import Any, Dict, List
import logging
import atexit
from constants import LOGGING_DIR
import os

class ParquetRequestLogger:
    """
    Thread-safe request logger that writes to Parquet files.
    Automatically handles file rotation and cleanup.
    """
    def __init__(
        self,
        log_dir: str = os.path.join(LOGGING_DIR, "requests"),
        schema: pa.Schema = None,
        max_rows_per_file: int = 100000,
        # retention_days: int = 30
    ):
        self.log_dir = Path(log_dir)
        self.max_rows = max_rows_per_file
        # self.retention_days = retention_days
        self.lock = threading.Lock()
        self.buffered_rows: List[Dict[str, Any]] = []
        
        # Default schema if none provided
        self.schema = schema or pa.schema([
            ('timestamp', pa.timestamp('us')),
            ('method', pa.string()),
            ('path', pa.string()),
            ('status_code', pa.int32()),
            ('response_time_ms', pa.float32()),
            ('client_ip', pa.string())
        ])
        
        # Ensure log directory exists
        self.log_dir.mkdir(parents=True, exist_ok=True)
        
        # Register cleanup on exit
        atexit.register(self.flush)
        # atexit.register(self.cleanup_old_logs)
    
    def _get_current_filename(self) -> Path:
        """Generate current log file path"""
        return self.log_dir / "requests.parquet"
    
    def log_request(self, data: Dict[str, Any]) -> None:
        """
        Thread-safely log request data to Parquet file.
        Buffers rows and writes in batches for efficiency.
        """
        if 'timestamp' not in data:
            data['timestamp'] = datetime.utcnow()
            
        with self.lock:
            self.buffered_rows.append(data)
            
            if len(self.buffered_rows) >= self.max_rows:
                self.flush()
    
    def flush(self) -> None:
        """Write buffered rows to Parquet file"""
        if not self.buffered_rows:
            return
            
        try:
            # Convert buffered rows to Arrow Table
            table = pa.Table.from_pylist(self.buffered_rows, schema=self.schema)
            
            output_file = self._get_current_filename()
            
            # Append to existing file if it exists
            if output_file.exists():
                existing_table = pq.read_table(output_file)
                table = pa.concat_tables([existing_table, table])
            
            pq.write_table(
                table,
                output_file,
                compression='snappy',
                row_group_size=10000
            )
            
            # Clear buffer
            self.buffered_rows = []
            
        except Exception as e:
            logging.error(f"Failed to write Parquet file: {e}")
            # Keep data in buffer to retry on next flush
    
    # def cleanup_old_logs(self) -> None:
    #     """Remove log files older than retention_days"""
    #     retention_threshold = datetime.now() - timedelta(days=self.retention_days)
        
    #     for log_file in self.log_dir.glob("requests_*.parquet"):
    #         # Extract date from filename
    #         try:
    #             file_date = datetime.strptime(
    #                 log_file.stem.split('_')[1],
    #                 "%Y%m%d"
    #             )
    #             if file_date < retention_threshold:
    #                 log_file.unlink()
    #         except Exception as e:
    #             logging.error(f"Failed to process/delete old log file {log_file}: {e}")

    # def read_logs(self, start_date: datetime = None, end_date: datetime = None) -> pa.Table:
    #     """
    #     Read logs within date range from the single Parquet file.
    #     Returns an Arrow Table for efficient processing.
    #     """
    #     log_file = self._get_current_filename()
        
    #     if not log_file.exists():
    #         return pa.Table.from_pylist([], schema=self.schema)
            
    #     table = pq.read_table(log_file, columns=self.schema.names)
        
    #     if start_date or end_date:
    #         timestamps = table['timestamp']
    #         mask = pa.array([True] * len(timestamps))
            
    #         if start_date:
    #             mask = mask & (timestamps >= pa.scalar(start_date))
    #         if end_date:
    #             mask = mask & (timestamps <= pa.scalar(end_date))
                
    #         table = table.filter(mask)
            
    #     return table

# Example usage
# if __name__ == "__main__":
#     # Create logger with custom schema (optional)
#     custom_schema = pa.schema([
#         ('timestamp', pa.timestamp('us')),
#         ('method', pa.string()),
#         ('path', pa.string()),
#         ('status_code', pa.int32()),
#         ('response_time_ms', pa.float32()),
#         ('client_ip', pa.string()),
#         # Add any additional fields you need
#         ('user_id', pa.string()),
#         ('payload_size', pa.int64())
#     ])
    
#     logger = ParquetRequestLogger(schema=custom_schema)
    
#     # Log a sample request
#     request_data = {
#         "method": "GET",
#         "path": "/api/users",
#         "status_code": 200,
#         "response_time_ms": 45.5,
#         "client_ip": "192.168.1.100",
#         "user_id": "user123",
#         "payload_size": 1024
#     }
#     logger.log_request(request_data)
    
#     # Force flush to disk (normally happens automatically at buffer limit)
#     logger.flush()
    
#     # Read logs for analysis
#     last_week = datetime.now() - timedelta(days=7)
#     recent_logs = logger.read_logs(start_date=last_week)
#     print(f"Total records: {recent_logs.num_rows}")