package edu.uci.ics.genomix.type;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.lang.reflect.ParameterizedType;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.io.WritableComparable;

public class ExternalableTreeSet<T extends WritableComparable<T>> implements Writable {

    static FileManager manager;

    public static void setupManager(Configuration conf, Path workPath) throws IOException {
        manager = new FileManager(conf, workPath);
    }

    static private int countLimit = Integer.MAX_VALUE;

    public static void setCountLimit(int count) {
        countLimit = count;
    }

    private TreeSet<T> inMemorySet;
    private Path path;
    private boolean isChanged;

    public ExternalableTreeSet() {
        this(null);
    }

    protected ExternalableTreeSet(Path path) {
        inMemorySet = new TreeSet<T>();
        this.path = path;
        isChanged = false;
    }

    public boolean add(T t) {
        boolean contains = inMemorySet.add(t);
        if (contains) {
            isChanged = contains;
        }
        return contains;
    }

    public boolean remove(T t) {
        boolean contains = inMemorySet.remove(t);
        if (contains) {
            isChanged = contains;
        }
        return contains;
    }

    public SortedSet<T> rangeSearch(T lowKey, T highKey) {
        SortedSet<T> set = inMemorySet.subSet(lowKey, highKey);
        return set;
    }

    /**
     * Union with setB, make sure setB have already loaded before use.
     * 
     * @return true if the set changed.
     * @param setB
     */
    public boolean union(ExternalableTreeSet<T> setB) {
        boolean changed = inMemorySet.addAll(setB.inMemorySet);
        if (changed) {
            isChanged = true;
        }
        return changed;
    }

    protected static TreeSet<?> load(Path path) throws IOException, ClassNotFoundException {
        if (path == null) {
            return null;
        }
        InputStream fis = manager.getInputStream(path);
        ObjectInputStream ois = new ObjectInputStream(fis);
        TreeSet<?> set = (TreeSet<?>) ois.readObject();
        ois.close();
        return set;
    }

    protected static void save(Path path, final TreeSet<?> set) throws IOException {
        if (path == null) {
            return;
        }
        OutputStream fos = manager.getOutputStream(path);
        ObjectOutputStream oos = new ObjectOutputStream(fos);
        oos.writeObject(set);
        oos.close();
    }

    @SuppressWarnings("unchecked")
    @Override
    public void readFields(DataInput in) throws IOException {
        int size = in.readInt();
        inMemorySet.clear();
        if (size < countLimit) {
            path = null;
            for (int i = 0; i < size; ++i) {
                Class<T> clazz = (Class<T>) ((ParameterizedType) getClass().getGenericSuperclass())
                        .getActualTypeArguments()[0];;
                T t;
                try {
                    t = clazz.newInstance();
                } catch (InstantiationException e) {
                    throw new IOException(e);
                } catch (IllegalAccessException e) {
                    throw new IOException(e);
                }
                t.readFields(in);
                inMemorySet.add(t);
            }
        } else {
            if (path == null) {
                path = new Path(in.readUTF());
            }
            try {
                inMemorySet = (TreeSet<T>) load(path);
            } catch (ClassNotFoundException e) {
                throw new IOException(e);
            }
        }
        isChanged = false;
    }

    @Override
    public void write(DataOutput out) throws IOException {
        out.writeInt(inMemorySet.size());
        if (inMemorySet.size() < countLimit) {
            for (T t : inMemorySet) {
                t.write(out);
            }
        } else {
            if (path == null) {
                path = manager.createFile();
            }
            if (isChanged) {
                save(path, inMemorySet);
            }
            out.writeUTF(path.toString());
        }
        isChanged = false;
    }

    public void destroy() throws IOException {
        manager.deleteFile(path);
    }

    protected static class FileManager {
        private FileSystem fs;
        private Path workPath;
        private HashMap<Path, OutputStream> allocatedPath;

        public FileManager(Configuration conf, Path workingPath) throws IOException {
            fs = FileSystem.get(conf);
            workPath = workingPath;
            allocatedPath = new HashMap<Path, OutputStream>();
        }

        protected final SimpleDateFormat simpleDateFormat = new SimpleDateFormat("ddMMyy-hhmmssSS");

        public synchronized Path createFile() throws IOException {
            Path path;
            do {
                path = new Path(workPath, simpleDateFormat.format(new Date()));
            } while (fs.exists(path));
            allocatedPath.put(path, fs.create(path, (short) 1));
            return path;
        }

        public synchronized void deleteFile(Path path) throws IOException {
            if (path != null && allocatedPath.containsKey(path)) {
                fs.delete(path, true);
                allocatedPath.remove(path);
            }
        }

        public OutputStream getOutputStream(Path path) throws IOException {
            if (!allocatedPath.containsKey(path)) {
                throw new IOException("File not exist:" + path);
            }
            return allocatedPath.get(path);
        }

        public InputStream getInputStream(Path path) throws IOException {
            if (!allocatedPath.containsKey(path)) {
                throw new IOException("File not exist:" + path);
            }
            return fs.open(path);
        }
    }
}
